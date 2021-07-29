# Standard imports
from glob import glob
import json

# Third-party imports
from netCDF4 import Dataset
import numpy as np

class Input:
    """Extracts and stores reach-level FLPE algorithm data.
    
    Attributes
    ----------
    alg_dict: dict
        dictionary of algorithm data stored by algorithm name as numpy arrays
    alg_dir: Path
        path to reach-level FLPE algorithm data
    basin_dict: dict
        dict of reach_ids and SoS file needed to process entire basin of data
    sos_dict: dict
        dictionary of SoS data
    sos_dir: Path
        path to SoS data    

    Methods
    -------
    extract_alg()
        extracts and stores reach-level FLPE algorithm data
    extract_sos()
        extracts and stores SoS data
    __get_ids(self, basin_json):
        Extract reach identifiers and store in basin_dict
    """

    def __init__(self, alg_dir, basin_json, sos_dir):
        """
        Parameters
        ----------
        alg_dir: Path
            path to reach-level FLPE algorithm data
        basin_json: Path
            path to basin JSON file
        sos_dir: Path
            path to SoS data
        """

        self.alg_dict = {
            "geobam": {},
            "hivdi": {},
            "metroman": {},
            "momma": {},
            "sad": {}
        }
        self.basin_dict = self.__get_ids(basin_json)
        self.alg_dir = alg_dir
        self.sos_dict = {}
        self.sos_dir = sos_dir

    def extract_sos(self):
        """Extracts and stores SoS data in sos_dict.

        TODO: Implement
        
        Parameters
        ----------
        """

        raise NotImplementedError

    def extract_alg(self):
        """Extracts and stores reach-level FLPE algorithm data in alg_dict.

        TODO: 
        - Modify to store basin-level data instead of reach.
        - Locate HiVDI n parameter.
        - Locate MOMMA A0 parameter.
        """

        reach_id = self.basin_dict["reach_ids"][0]    ## TODO

        # geobam
        gb_file = self.alg_dir / "geobam" / f"{reach_id}_geobam.nc"
        gb = Dataset(gb_file, 'r', format="NETCDF4")
        self.alg_dict["geobam"]["q"] = self.__get_gb_data(gb, "logQ", True)
        self.alg_dict["geobam"]["n"] = self.__get_gb_data(gb, "logn_man", True)
        self.alg_dict["geobam"]["a0"] = self.__get_gb_data(gb, "A0", False)
        gb.close()

        # hivdi
        hivdi_file = self.alg_dir / "hivdi" / f"{reach_id}_hivdi.nc"
        hv = Dataset(hivdi_file, 'r', format="NETCDF4")
        self.alg_dict["hivdi"]["q"] = hv["reach"]["Q"][:].filled(np.nan)
        self.alg_dict["hivdi"]["n"] = hv["reach"]["alpha"][:].filled(np.nan)  ## TODO
        self.alg_dict["hivdi"]["a0"] = hv["reach"]["A0"][:].filled(np.nan)
        hv.close()

        # momma
        momma_file = self.alg_dir / "momma" / f"{reach_id}_momma.nc"
        mo = Dataset(momma_file, 'r', format="NETCDF4")
        self.alg_dict["momma"]["q"] = mo["Q"][:].filled(np.nan)
        self.alg_dict["momma"]["n"] = mo["n"][:].filled(np.nan)
        self.alg_dict["momma"]["a0"] = mo["Y"][:].filled(np.nan)    ## TODO
        mo.close()

        # sad
        sad_file = self.alg_dir / "sad" / f"{reach_id}_sad.nc"
        sd = Dataset(sad_file, 'r', format="NETCDF4")
        self.alg_dict["sad"]["q"] = sd["Qa"][:].filled(np.nan)
        self.alg_dict["sad"]["n"] = sd["n"][:].filled(np.nan)
        self.alg_dict["sad"]["a0"] = sd["A0"][:].filled(np.nan)    ## TODO
        sd.close()

        # metroman
        metroman_file = [ name for name in glob(str(self.alg_dir / "metroman" / f"*{reach_id}*_metroman.nc"))][0]
        mm = Dataset(metroman_file, 'r', format="NETCDF4")
        index = np.where(mm["reach_id"][:] == reach_id)
        self.alg_dict["metroman"]["q"] = mm["allq"][index].filled(np.nan)
        self.alg_dict["metroman"]["n"] = mm["nahat"][index].filled(np.nan)
        self.alg_dict["metroman"]["a0"] = mm["A0hat"][index].filled(np.nan)
        mm.close()        

    def __get_gb_data(self, gb, group, logged):
        """Return geoBAM data as a numpy array.
        
        Parameters
        ----------
        gb: netCDF4.Dataset
            NetCDF file dataset to extract discharge time series
        group: str
            string name of group to access chains
        logged: bool
            boolean indicating if result is logged
        """

        chain1 = gb[group]["mean_chain1"][:].filled(np.nan)
        chain2 = gb[group]["mean_chain2"][:].filled(np.nan)
        chain3 = gb[group]["mean_chain3"][:].filled(np.nan)
        chains = np.vstack((chain1, chain2, chain3))
        
        if logged:
            return np.exp(np.nanmean(chains, axis=0))
        else:
            return np.nanmean(chains, axis=0)

    def __get_ids(self, basin_json):
        """Extract reach identifiers and return dictionary.
        
        Dictionary is organized with a key of reach identifier and a value of
        SoS file as a Path object.

        TODO: Organize with key of basin number and value of dict with reach ids
        and sos file.
        """

        # index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
        index = 13
        with open(basin_json) as json_file:
            data = json.load(json_file)

        return {
            "basin_id" : int(str(data[index]["reach_id"])[1:6]),
            "reach_ids" : [data[index]["reach_id"]],    ## TODO
            "sos" : data[index]["sos"]
        }