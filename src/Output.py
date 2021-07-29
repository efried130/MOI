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

    def __init__(self, basin_dict, out_dir, integ_dict):
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

    def write_output(self):
        """Write data stored to NetCDF file labelled with basin id.
        
        TODO: 
        - Add optional attribute metadata like valid
        - Decide on how to store output of FLPE integration (how to identify
        which data corresponds to which algorithm)
        """

        # NetCDF file creation
        out_file = self.out_dir / f"{str(self.basin_dict['basin_id'])}_integrator.nc"
        out = Dataset(out_file, 'w', format="NETCDF4")
        out.production_date = datetime.now().strftime('%d-%b-%Y %H:%M:%S')
        
        # Dimensions and coordinate variables
        out.createDimension("nr", len(self.basin_dict["reach_ids"]))
        out.createDimension("nflpe", len(self.stage_estimate["flpe"].keys()))
        nr = out.createVariable("nr", "i4", ("nr",))
        nr.units = "reaches"
        nr.long_name = "number of reaches"
        nr[:] = range(1, len(self.basin_dict["reach_ids"]) + 1)
        # nf = out.createVariable("nflpe", "i4", ("nflpe",))
        # nf.units = "algorithms"
        # nf.long_name = "number of reach-level FLPE algorithms"
        # nf[:] = range(1, len(self.stage_estimate["flpe"].keys()) + 1)    ## TODO

        # Data
        # reach identifiers
        reach_id_v = out.createVariable("reach_id", "i8")
        reach_id_v.long_name = "reach ID from prior river database"
        reach_id_v.comment = "Unique reach identifier from the prior river " \
            + "database. The format of the identifier is CBBBBBRRRRT, where " \
            + "C=continent, B=basin, R=reach, T=type."
        reach_id_v[:] = np.array(self.basin_dict["reach_ids"], dtype=int)
        
        # mean discharge
        q_var = out.createVariable("Qmean", "f4", ("nr",),
            fill_value=self.FILL_VALUE)
        q_var[:] = self.stage_estimate["q_mean"][:]

        # FLPE results
        gb_var = out.createVariable("geobam", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        gb_var[:] = self.stage_estimate["flpe"]["geobam"]

        hv_var = out.createVariable("hivdi", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        hv_var[:] = self.stage_estimate["flpe"]["hivdi"]

        mm_var = out.createVariable("metroman", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        mm_var[:] = self.stage_estimate["flpe"]["metroman"]

        mo_var = out.createVariable("momma", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        mo_var[:] = self.stage_estimate["flpe"]["momma"]

        sd_var = out.createVariable("sad", "f4", ("nr"),
            fill_value=self.FILL_VALUE)
        sd_var[:] = self.stage_estimate["flpe"]["sad"]

        out.close()