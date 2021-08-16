# Standard imports
from datetime import datetime

# Third-party imports
from netCDF4 import Dataset
import numpy as np

class Output:
    """Writes integration results stored in integ_dict to NetCDF file.
    
    Attributes
    ----------
    basin_dict: dict
        dict of reach_ids and SoS file needed to process entire basin of data
    FILL_VALUE: float
        Float fill value for missing data
    out_dir: Path
        path to output dir
    stage_estimate: dict
        dict of integrator estimate data

    Methods
    -------
    write_output()
        Write data stored to NetCDF file labelled with basin id
    """
    
    FILL_VALUE = -999999999999

    def __init__(self, basin_dict, out_dir, integ_dict, alg_dict, obs_dict):
        """
        Parameters
        ----------
        basin_dict: dict
            dict of reach_ids and SoS file needed to process entire basin of data
        out_dir: Path
            path to output dir
        integ_dict: dict
            dict of integrator estimate data
        """

        self.basin_dict = basin_dict
        self.out_dir = out_dir
        self.stage_estimate = integ_dict
        self.alg_dict = alg_dict
        self.obs_dict = obs_dict

    def write_output(self):
        """Write data stored to NetCDF files for each reach
        
        TODO: 
        - Add optional attribute metadata like valid
        - Should variable output be (nr, nt)?
        - Storage of FLPE algorithms
        - Storage of actual Integrator processing
        """

        fillvalue = -999999999999

        for reach in self.basin_dict['reach_ids']:

             iDelete=self.obs_dict[reach]['iDelete']
             shape_iDelete=np.shape(iDelete)
             nDelete=shape_iDelete[1]
             iInsert=iDelete-np.arange(nDelete)
             iInsert=np.reshape(iInsert,[nDelete,]) 

             self.obs_dict[reach]['nt'] += nDelete

             #print(self.alg_dict['geobam'][reach]['integrator']['q'])
             #print(self.alg_dict['geobam'][reach]['integrator']['q'].shape)

             self.alg_dict['geobam'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['geobam'][reach]['integrator']['q'],iInsert,fillvalue,1)

             # NetCDF file creation
             out_file = self.out_dir / f"{reach}_integrator.nc"
             out = Dataset(out_file, 'w', format="NETCDF4")
             out.production_date = datetime.now().strftime('%d-%b-%Y %H:%M:%S')

             # Dimensions and coordinate variables
             out.createDimension("nt", self.obs_dict[reach]['nt'] )
             nt = out.createVariable("nt", "i4", ("nt",))
             nt.units = "time steps"
             nt[:] = range(self.obs_dict[reach]['nt'])

             # Geobam
             gb = out.createGroup("geobam")
             gbq  = out.createVariable("geobam/q", "f8", ("nt",), fill_value=fillvalue)
             gbq[:] = self.alg_dict['geobam'][reach]['integrator']['q']
             gb_a0  = out.createVariable("geobam/a0", "f8", fill_value=fillvalue)
             gb_a0[:] = self.alg_dict['geobam'][reach]['integrator']['a0']
             gb_n  = out.createVariable("geobam/n", "f8", fill_value=fillvalue)
             gb_n[:] = self.alg_dict['geobam'][reach]['integrator']['n']
             gb_qbar_stage1  = out.createVariable("geobam/qbar_reachScale", "f8", fill_value=fillvalue)
             gb_qbar_stage1[:] = self.alg_dict['geobam'][reach]['qbar']
             gb_qbar_stage2  = out.createVariable("geobam/qbar_basinScale", "f8", fill_value=fillvalue)
             gb_qbar_stage2[:] = self.alg_dict['geobam'][reach]['integrator']['qbar']

             out.close()

        '''
        

        # out.createDimension("nflpe", len(self.stage_estimate["flpe"].keys()))
        # nf = out.createVariable("nflpe", "i4", ("nflpe",))
        # nf.units = "algorithms"
        # nf.long_name = "number of reach-level FLPE algorithms"
        # nf[:] = range(1, len(self.stage_estimate["flpe"].keys()) + 1)    ## TODO Figure out time series and FLPE storage

        # Data
        # reach identifiers
        reach_id_v = out.createVariable("reach_id", "i8", ("nr",))
        reach_id_v.long_name = "reach ID from prior river database"
        reach_id_v.comment = "Unique reach identifier from the prior river " \
            + "database. The format of the identifier is CBBBBBRRRRT, where " \
            + "C=continent, B=basin, R=reach, T=type."
        reach_id_v[:] = np.array(self.basin_dict["reach_ids"], dtype=int)

        # pre integrator mean discharge
        pre_q_var = out.createVariable("pre_Qmean", "f4", ("nr",),
            fill_value=self.FILL_VALUE)
        pre_q_var[:] = np.random.uniform(size=len(nr[:]))

        # mean discharge
        q_var = out.createVariable("Qmean", "f4", ("nr",),
            fill_value=self.FILL_VALUE)
        # q_var[:] = self.stage_estimate["q_mean"][:]    ## TODO Store processing
        q_var[:] = np.random.uniform(size=len(nr[:]))

        # FLPE results
        gb_var = out.createVariable("geobam", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        # gb_var[:] = self.stage_estimate["flpe"]["geobam"]    ## TODO Store processing
        gb_var[:] = np.random.uniform(size=len(nr[:]))

        hv_var = out.createVariable("hivdi", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        # hv_var[:] = self.stage_estimate["flpe"]["hivdi"]    ## TODO Store processing
        hv_var[:] = np.random.uniform(size=len(nr[:]))

        mm_var = out.createVariable("metroman", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        # mm_var[:] = self.stage_estimate["flpe"]["metroman"]    ## TODO Store processing
        mm_var[:] = np.random.uniform(size=len(nr[:]))

        mo_var = out.createVariable("momma", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        # mo_var[:] = self.stage_estimate["flpe"]["momma"]    ## TODO Store processing
        mo_var[:] = np.random.uniform(size=len(nr[:]))

        sd_var = out.createVariable("sad", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        # sd_var[:] = self.stage_estimate["flpe"]["sad"]    ## TODO Store processing
        sd_var[:] = np.random.uniform(size=len(nr[:]))

        '''
