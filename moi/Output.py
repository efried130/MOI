# Standard imports
from datetime import datetime

# Third-party imports
from netCDF4 import Dataset
import numpy as np
import shutil

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

    def __init__(self, basin_dict, out_dir, integ_dict, alg_dict, obs_dict, sword_dir):
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
        self.sword_dir = sword_dir

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

             self.alg_dict['geobam'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['geobam'][reach]['integrator']['q'],iInsert,fillvalue,1)

             self.alg_dict['hivdi'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['hivdi'][reach]['integrator']['q'],iInsert,fillvalue,1)

             self.alg_dict['metroman'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['metroman'][reach]['integrator']['q'],iInsert,fillvalue,1)

             self.alg_dict['momma'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['momma'][reach]['integrator']['q'],iInsert,fillvalue,1)

             self.alg_dict['sad'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['sad'][reach]['integrator']['q'],iInsert,fillvalue,1)

             self.alg_dict['sic4dvar'][reach]['integrator']['q']=np.insert( \
                   self.alg_dict['sic4dvar'][reach]['integrator']['q'],iInsert,fillvalue,1)

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
             gbq[:] = np.nan_to_num(self.alg_dict['geobam'][reach]['integrator']['q'], copy=True, nan=fillvalue)
             
             gb_a0  = out.createVariable("geobam/a0", "f8", fill_value=fillvalue)
             gb_a0[:] = np.nan_to_num(self.alg_dict['geobam'][reach]['integrator']['a0'], copy=True, nan=fillvalue)
             
             gb_n  = out.createVariable("geobam/n", "f8", fill_value=fillvalue)
             gb_n[:] = np.nan_to_num(self.alg_dict['geobam'][reach]['integrator']['n'], copy=True, nan=fillvalue)
             
             gb_qbar_stage1  = out.createVariable("geobam/qbar_reachScale", "f8", fill_value=fillvalue)
             gb_qbar_stage1[:] = np.nan_to_num(self.alg_dict['geobam'][reach]['qbar'], copy=True, nan=fillvalue)
             
             gb_qbar_stage2  = out.createVariable("geobam/qbar_basinScale", "f8", fill_value=fillvalue)
             gb_qbar_stage2[:] = np.nan_to_num(self.alg_dict['geobam'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             # hivdi
             hv = out.createGroup("hivdi")
             hvq  = out.createVariable("hivdi/q", "f8", ("nt",), fill_value=fillvalue)
             hvq[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['integrator']['q'], copy=True, nan=fillvalue)

             hv_Abar = out.createVariable("hivdi/Abar", "f8", fill_value=fillvalue)
             hv_Abar[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['integrator']['Abar'], copy=True, nan=fillvalue)

             hv_alpha = out.createVariable("hivdi/alpha", "f8", fill_value=fillvalue)
             hv_alpha[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['integrator']['alpha'], copy=True, nan=fillvalue)

             hv_beta = out.createVariable("hivdi/beta", "f8", fill_value=fillvalue)
             hv_beta[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['integrator']['beta'], copy=True, nan=fillvalue)

             hv_qbar_stage1  = out.createVariable("hivdi/qbar_reachScale", "f8", fill_value=fillvalue)
             hv_qbar_stage1[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['qbar'], copy=True, nan=fillvalue)
             
             hv_qbar_stage2  = out.createVariable("hivdi/qbar_basinScale", "f8", fill_value=fillvalue)
             hv_qbar_stage2[:] = np.nan_to_num(self.alg_dict['hivdi'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             # metroman
             mm = out.createGroup("metroman")
             mmq  = out.createVariable("metroman/q", "f8", ("nt",), fill_value=fillvalue)
             mmq[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['q'], copy=True, nan=fillvalue)

             mm_Abar = out.createVariable("metroman/Abar", "f8", fill_value=fillvalue)
             mm_Abar[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['a0'], copy=True, nan=fillvalue)

             mm_na = out.createVariable("metroman/na", "f8", fill_value=fillvalue)
             mm_na[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['na'], copy=True, nan=fillvalue)

             mm_x1 = out.createVariable("metroman/x1", "f8", fill_value=fillvalue)
             mm_x1[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['x1'], copy=True, nan=fillvalue)

             mm_qbar_stage1  = out.createVariable("metroman/qbar_reachScale", "f8", fill_value=fillvalue)
             mm_qbar_stage1[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['qbar'], copy=True, nan=fillvalue)
             
             mm_qbar_stage2  = out.createVariable("metroman/qbar_basinScale", "f8", fill_value=fillvalue)
             mm_qbar_stage2[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             mm_sbQ_rel = out.createVariable("metroman/sbQ_rel", "f8", fill_value=fillvalue)
             mm_sbQ_rel[:] = np.nan_to_num(self.alg_dict['metroman'][reach]['integrator']['sbQrel'], copy=True, nan=fillvalue)

             # momma
             mo = out.createGroup("momma")
             moq  = out.createVariable("momma/q", "f8", ("nt",), fill_value=fillvalue)
             moq[:] = np.nan_to_num(self.alg_dict['momma'][reach]['integrator']['q'], copy=True, nan=fillvalue)

             mo_B = out.createVariable("momma/B", "f8", fill_value=fillvalue)
             mo_B[:] = np.nan_to_num(self.alg_dict['momma'][reach]['integrator']['B'], copy=True, nan=fillvalue)

             mo_H = out.createVariable("momma/H", "f8", fill_value=fillvalue)
             mo_H[:] = np.nan_to_num(self.alg_dict['momma'][reach]['integrator']['H'], copy=True, nan=fillvalue)

             mo_Save = out.createVariable("momma/Save", "f8", fill_value=fillvalue)
             mo_Save[:] = np.nan_to_num(self.alg_dict['momma'][reach]['integrator']['Save'], copy=True, nan=fillvalue)

             mo_qbar_stage1  = out.createVariable("momma/qbar_reachScale", "f8", fill_value=fillvalue)
             mo_qbar_stage1[:] = np.nan_to_num(self.alg_dict['momma'][reach]['qbar'], copy=True, nan=fillvalue)
             
             mo_qbar_stage2  = out.createVariable("momma/qbar_basinScale", "f8", fill_value=fillvalue)
             mo_qbar_stage2[:] = np.nan_to_num(self.alg_dict['momma'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             #sad
             sad=out.createGroup("sad")
             sadq = out.createVariable("sad/q", "f8", ("nt",), fill_value=fillvalue)
             sadq[:] = np.nan_to_num(self.alg_dict['sad'][reach]['integrator']['q'], copy=True, nan=fillvalue)

             sad_n = out.createVariable("sad/n", "f8", fill_value=fillvalue)
             sad_n[:] = np.nan_to_num(self.alg_dict['sad'][reach]['integrator']['n'], copy=True, nan=fillvalue)

             sad_a0 = out.createVariable("sad/a0", "f8", fill_value=fillvalue)
             sad_a0[:] = np.nan_to_num(self.alg_dict['sad'][reach]['integrator']['a0'], copy=True, nan=fillvalue)

             sad_qbar_stage1  = out.createVariable("sad/qbar_reachScale", "f8", fill_value=fillvalue)
             sad_qbar_stage1[:] = np.nan_to_num(self.alg_dict['sad'][reach]['qbar'], copy=True, nan=fillvalue)
             
             sad_qbar_stage2  = out.createVariable("sad/qbar_basinScale", "f8", fill_value=fillvalue)
             sad_qbar_stage2[:] = np.nan_to_num(self.alg_dict['sad'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             #sic4dvar
             sic4dvar=out.createGroup("sic4dvar")
             sic4dvarq = out.createVariable("sic4dvar/q", "f8", ("nt",), fill_value=fillvalue)
             sic4dvarq[:] = np.nan_to_num(self.alg_dict['sic4dvar'][reach]['integrator']['q'], copy=True, nan=fillvalue)

             sic4dvar_n = out.createVariable("sic4dvar/n", "f8", fill_value=fillvalue)
             sic4dvar_n[:] = np.nan_to_num(self.alg_dict['sic4dvar'][reach]['integrator']['n'], copy=True, nan=fillvalue)

             sic4dvar_a0 = out.createVariable("sic4dvar/a0", "f8", fill_value=fillvalue)
             sic4dvar_a0[:] = np.nan_to_num(self.alg_dict['sic4dvar'][reach]['integrator']['a0'], copy=True, nan=fillvalue)

             sic4dvar_qbar_stage1  = out.createVariable("sic4dvar/qbar_reachScale", "f8", fill_value=fillvalue)
             sic4dvar_qbar_stage1[:] = np.nan_to_num(self.alg_dict['sic4dvar'][reach]['qbar'], copy=True, nan=fillvalue)
             
             sic4dvar_qbar_stage2  = out.createVariable("sic4dvar/qbar_basinScale", "f8", fill_value=fillvalue)
             sic4dvar_qbar_stage2[:] = np.nan_to_num(self.alg_dict['sic4dvar'][reach]['integrator']['qbar'], copy=True, nan=fillvalue)

             out.close()

    def write_sword_output(self):
        """Make a new copy of the SWORD file, and write the Confluence estimates of the FLPs into the file.
           by Mike, September 2022
           """

        sword_src_file=self.sword_dir.joinpath(self.basin_dict['sword'])
        sword_dest_file=self.out_dir.joinpath(self.basin_dict['sword'])
        shutil.copy(sword_src_file,sword_dest_file)

        sword_dataset = Dataset(sword_dest_file,'a')

        reaches = sword_dataset['reaches']['reach_id'][:]

        for reach in self.basin_dict['reach_ids']:
            reach_ind = np.where(reaches == np.int(reach))
            #print(reach)
            #print(reach_ind)
 
            #1 bam 
            sword_dataset['reaches']['discharge_models']['unconstrained']['BAM']['Abar'][reach_ind]= \
                self.alg_dict['geobam'][reach]['integrator']['a0']
            sword_dataset['reaches']['discharge_models']['unconstrained']['BAM']['n'][reach_ind]= \
                self.alg_dict['geobam'][reach]['integrator']['n']
            #2 hivdi
            sword_dataset['reaches']['discharge_models']['unconstrained']['HiVDI']['Abar'][reach_ind]=\
                self.alg_dict['hivdi'][reach]['integrator']['Abar']
            sword_dataset['reaches']['discharge_models']['unconstrained']['HiVDI']['alpha'][reach_ind]=\
                self.alg_dict['hivdi'][reach]['integrator']['alpha']
            sword_dataset['reaches']['discharge_models']['unconstrained']['HiVDI']['beta'][reach_ind]=\
                self.alg_dict['hivdi'][reach]['integrator']['beta']
            #3 metroman
            sword_dataset['reaches']['discharge_models']['unconstrained']['MetroMan']['Abar'][reach_ind]= \
                self.alg_dict['metroman'][reach]['integrator']['a0']
            sword_dataset['reaches']['discharge_models']['unconstrained']['MetroMan']['ninf'][reach_ind]= \
                self.alg_dict['metroman'][reach]['integrator']['na']
            sword_dataset['reaches']['discharge_models']['unconstrained']['MetroMan']['p'][reach_ind]= \
                self.alg_dict['metroman'][reach]['integrator']['x1']
            #4 momma
            sword_dataset['reaches']['discharge_models']['unconstrained']['MOMMA']['B'][reach_ind]= \
                self.alg_dict['momma'][reach]['integrator']['B']
            sword_dataset['reaches']['discharge_models']['unconstrained']['MOMMA']['H'][reach_ind]= \
                self.alg_dict['momma'][reach]['integrator']['H']
            sword_dataset['reaches']['discharge_models']['unconstrained']['MOMMA']['Save'][reach_ind]= \
                self.alg_dict['momma'][reach]['integrator']['Save']
            #5 sads
            sword_dataset['reaches']['discharge_models']['unconstrained']['SADS']['Abar'][reach_ind]= \
                self.alg_dict['sad'][reach]['integrator']['a0']
            sword_dataset['reaches']['discharge_models']['unconstrained']['SADS']['n'][reach_ind]= \
                self.alg_dict['sad'][reach]['integrator']['n']
            #6 sicvdvar - note this is not included in SWORD v11 - so can't add these, at the moment
 

        sword_dataset.close()
