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

    def __init__(self, alg_dir, basin_json, sos_dir, swot_dir):
        """
        Parameters
        ----------
        alg_dir: Path
            path to reach-level FLPE algorithm data
        basin_json: Path
            path to basin JSON file
        sos_dir: Path
            path to SoS data
        swot_dir: Path
            path to SWOT data
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
        self.swot_dir = swot_dir

    def extract_sos(self):
        """Extracts and stores SoS data in sos_dict.

        TODO: Implement
        
        Parameters
        ----------
        """

        raise NotImplementedError

    def extract_swot(self):
 
        self.obs_dict={}
#        self.obs_dict['h']={}
#        self.obs_dict['w']={}
#        self.obs_dict['S']={}
#        self.obs_dict['dA']={}

        for reach in self.basin_dict['reach_ids']:
             swotfile=self.swot_dir.joinpath(reach+'_SWOT.nc')
             swot_dataset = Dataset(swotfile)

             self.obs_dict[reach]={}
             nt=len(swot_dataset.dimensions['nt'])
             self.obs_dict[reach]['nt']=nt
             self.obs_dict[reach]['h']=swot_dataset["reach/wse"][0:nt].filled(np.nan)
             self.obs_dict[reach]['w']=swot_dataset["reach/width"][0:nt].filled(np.nan)
             self.obs_dict[reach]['S']=swot_dataset["reach/slope2"][0:nt].filled(np.nan)
             self.obs_dict[reach]['dA']=swot_dataset["reach/d_x_area"][0:nt].filled(np.nan)

             swot_dataset.close()

             #select observations that are NOT equal to the fill value
             iDelete=np.where(np.isnan(self.obs_dict[reach]['h']))
             shape_iDelete=np.shape(iDelete)
             nDelete=shape_iDelete[1]
             self.obs_dict[reach]['h']=np.delete(self.obs_dict[reach]['h'],iDelete,0)
             self.obs_dict[reach]['w']=np.delete(self.obs_dict[reach]['w'],iDelete,0)
             self.obs_dict[reach]['S']=np.delete(self.obs_dict[reach]['S'],iDelete,0)
             self.obs_dict[reach]['dA']=np.delete(self.obs_dict[reach]['dA'],iDelete,0)

             self.obs_dict[reach]['nt'] -= nDelete

    def extract_alg(self):
        """Extracts and stores reach-level FLPE algorithm data in alg_dict.

        TODO: 
        - Consider storage of basin-level data in a DataFrame if possible.
        - Locate HiVDI n parameter.
        - Locate MOMMA A0 parameter.
        - Add in MetroMan results once available.
        - Add in SAD results once available.
        """

        reach_ids = self.basin_dict["reach_ids"]
        for r_id in reach_ids:
            gb_file = self.alg_dir / "geobam" / f"{r_id}_geobam.nc"
            hv_file = self.alg_dir / "hivdi" / f"{r_id}_hivdi.nc"
            mo_file = self.alg_dir / "momma" / f"{r_id}_momma.nc"
            sd_file = self.alg_dir / "sad" / f"{r_id}_sad.nc"      ## TODO wait on SAD results
            mm_file = glob(str(self.alg_dir / "metroman" / f"*{r_id}*_metroman.nc"))    ## TODO hold until results are available

            if gb_file.exists() and hv_file.exists() and mo_file.exists():    ## TODO add in SAD and MetroMan files

                self.__extract_valid(r_id, gb_file, hv_file, mo_file, sd_file, mm_file)
            else:
                self.__indicate_no_data(r_id)

    def __extract_valid(self, r_id, gb_file, hv_file, mo_file, sd_file, mm_file):
        """ Extract valid data from the output of each reach-level FLPE alg.

        TODO: Metroman results

        Parameters
        ----------
        r_id: str
            Unique reach identifier
        gb_file: Path
            Path to geoBAM results file
        hv_file: Path
            Path to HiVDI results file
        mo_file: Path
            Path to MOMMA results file
        sd_file: Path
            Path to SAD results file
        mm_file: Path
            Path to MetroMan results file
        """

        # geobam
        gb = Dataset(gb_file, 'r', format="NETCDF4")
        self.alg_dict["geobam"][r_id] = {
            "q": np.array(self.__get_gb_data(gb, "logQ", True)),
            "n": np.array(self.__get_gb_data(gb, "logn_man", True)),
            "a0": np.array(self.__get_gb_data(gb, "A0", False))
        }
        gb.close()

        # hivdi
        hv = Dataset(hv_file, 'r', format="NETCDF4")
        self.alg_dict["hivdi"][r_id] = {
            "q" : hv["reach"]["Q"][:].filled(np.nan),
            # "n" : hv["reach"]["alpha"][:].filled(np.nan),  ## TODO locate n
            "a0" : hv["reach"]["A0"][:].filled(np.nan)
        }
        hv.close()

        # momma
        mo = Dataset(mo_file, 'r', format="NETCDF4")
        self.alg_dict["momma"][r_id] = {
            "q" : mo["Q"][:].filled(np.nan),
            "n" : mo["n"][:].filled(np.nan)
            # "a0" : ...                                    ## TODO locate or calculate A0
        }
        mo.close()

        # sad
        # sd = Dataset(sd_file, 'r', format="NETCDF4")
        # self.alg_dict["sad"][r_id] = {
        #     "q" : sd["Qa"][:].filled(np.nan),
        #     "n" : sd["n"][:].filled(np.nan),
        #     "a0" : sd["A0"][:].filled(np.nan)
        # }
        # sd.close()

        # metroman    ## TODO hold until results are available
        # mm = Dataset(mm_file, 'r', format="NETCDF4")
        # index = np.where(mm["reach_id"][:] == r_id)
        # self.alg_dict["metroman"][r_id] = {
        #     "q" : mm["allq"][index].filled(np.nan),
        #     "n" : mm["nahat"][index].filled(np.nan),
        #     "a0" : mm["A0hat"][index].filled(np.nan)
        # }
        # mm.close()

    def __indicate_no_data(self, r_id):
        """Indicate no data is available for the reach.

        TODO: Metroman results

        Parameters
        ----------
        r_id: str
            Unique reach identifier
        """

        self.alg_dict["geobam"][r_id] = {
            "q": np.nan,
            "n": np.nan,
            "a0": np.nan
        }

        # hivdi
        self.alg_dict["hivdi"][r_id] = {
            "q" : np.nan,
            # "n" : np.nan,  ## TODO locate n
            "a0" : np.nan
        }

        # momma
        self.alg_dict["momma"][r_id] = {
            "q" : np.nan,
            "n" : np.nan
            # "a0" : np.nan                                 ## TODO locate or calculate A0
        }

        # sad
        # self.alg_dict["sad"][r_id] = {
        #     "q" : np.nan,
        #     "n" : np.nan,
        #     "a0" : np.nan
        # }

        # MetroMan    ## TODO hold until results are available
        # self.alg_dict["metroman"][r_id] = {
        #     "q" : np.nan,
        #     "n" : np.nan,
        #     "a0" : np.nan
        # }

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
        """

        # index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
        index = 3
        with open(basin_json) as json_file:
            data = json.load(json_file)

        return {
            "basin_id" : int(data[index]["basin_id"]),    ## TODO
            "reach_ids" : data[index]["reach_id"],    ## TODO
            "sos" : data[index]["sos"],
            "sword": data[index]["sword"]
        }
