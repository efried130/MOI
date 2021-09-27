# Standard imports
from glob import glob
from pathlib import Path
import warnings

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

    def __init__(self, alg_dir, sos_dir, swot_dir, basin_data):
        """
        Parameters
        ----------
        alg_dir: Path
            path to reach-level FLPE algorithm data
        sos_dir: Path
            path to SoS data
        swot_dir: Path
            path to SWOT data
        basin_data: dict
            dict of reach_ids and SoS file needed to process entire basin of data
        """

        self.alg_dict = {
            "geobam": {},
            "hivdi": {},
            "metroman": {},
            "momma": {},
            "sad": {},
            "sic4dvar": {}
        }
        self.basin_dict = basin_data
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
             self.obs_dict[reach]['h']=np.delete(self.obs_dict[reach]['h'],iDelete,0)
             self.obs_dict[reach]['w']=np.delete(self.obs_dict[reach]['w'],iDelete,0)
             self.obs_dict[reach]['S']=np.delete(self.obs_dict[reach]['S'],iDelete,0)
             self.obs_dict[reach]['dA']=np.delete(self.obs_dict[reach]['dA'],iDelete,0)

             self.obs_dict[reach]['iDelete']=iDelete

             Smin=1.7e-5
             self.obs_dict[reach]['S'][self.obs_dict[reach]['S']<Smin]=\
                np.putmask(self.obs_dict[reach]['S'],self.obs_dict[reach]['S']<Smin,Smin)

             #Obs.S[Obs.S<Smin]=putmask(Obs.S,Obs.S<Smin,Smin) #limit slopes to a minimum value

             shape_iDelete=np.shape(iDelete)
             nDelete=shape_iDelete[1]
             self.obs_dict[reach]['nt'] -= nDelete

    def extract_alg(self):
        """Extracts and stores reach-level FLPE algorithm data in alg_dict."""

        reach_ids = self.basin_dict["reach_ids"]
        for r_id in reach_ids:

            gb_file = self.alg_dir / "geobam" / f"{r_id}_geobam.nc"
            hv_file = self.alg_dir / "hivdi" / f"{r_id}_hivdi.nc"
            mo_file = self.alg_dir / "momma" / f"{r_id}_momma.nc"
            sd_file = self.alg_dir / "sad" / f"{r_id}_sad.nc"
            sv_file = self.alg_dir / "sic4dvar" / f"{r_id}_sic4dvar.nc"
            mm_file = glob(str(self.alg_dir / "metroman" / f"*{r_id}*_metroman.nc"))    
            mm_file = Path(mm_file[0]) 

            if gb_file.exists() and hv_file.exists() and mm_file.exists() \
                and mo_file.exists() and sd_file.exists() and sv_file.exists():

                self.__extract_valid(r_id, gb_file, hv_file, mo_file, sd_file, mm_file, sv_file)
            else:
                self.__indicate_no_data(r_id)

    def __extract_valid(self, r_id, gb_file, hv_file, mo_file, sd_file, mm_file, sv_file):
        """ Extract valid data from the output of each reach-level FLPE alg.
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
        sv_file: Path
            Path to SIC4DVar results file
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
            "alpha" : hv["reach"]["alpha"][:].filled(np.nan),  
            "beta" : hv["reach"]["beta"][:].filled(np.nan),  
            "a0" : hv["reach"]["A0"][:].filled(np.nan)
        }
        hv.close()

        # momma
        mo = Dataset(mo_file, 'r', format="NETCDF4")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.alg_dict["momma"][r_id] = {
                "q" : mo["Q"][:].filled(np.nan),
                "B" : mo["zero_flow_stage"][:].filled(np.nan),
                "H" : mo["bankfull_stage"][:].filled(np.nan),                                  
                "Save" : np.nanmean(mo["slope"][:].filled(np.nan))
            }
        mo.close()

        # sad
        sd = Dataset(sd_file, 'r', format="NETCDF4")
        self.alg_dict["sad"][r_id] = {
            "q" : sd["Qa"][:].filled(np.nan),
            "n" : sd["n"][:].filled(np.nan),
            "a0" : sd["A0"][:].filled(np.nan)
        }
        sd.close()

        # metroman    
        mm = Dataset(mm_file, 'r', format="NETCDF4")
        index = np.where(mm["reach_id"][:] == int(r_id))
        self.alg_dict["metroman"][r_id] = {
             "q" : mm["allq"][index].filled(np.nan),
             "na" : mm["nahat"][index].filled(np.nan),
             "x1" : mm["x1hat"][index].filled(np.nan),
             "a0" : mm["A0hat"][index].filled(np.nan)
        }
        mm.close()

        # sic4dvar
        sv = Dataset(sv_file, 'r', format="NETCDF4")
        self.alg_dict["sic4dvar"][r_id] = {
            #"q31": sv["Qalgo31"][:].filled(np.nan),#unclear which of these to use
            "q": sv["Qalgo31"][:].filled(np.nan),
            "q5": sv["Qalgo5"][:].filled(np.nan),
            "n": sv["n"][:].filled(np.nan),
            "a0": sv["A0"][:].filled(np.nan)
        }
        sv.close()

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
            "alpha" : np.nan,  
            "beta" : np.nan,  
            "a0" : np.nan
        }

        # momma
        self.alg_dict["momma"][r_id] = {
            "q" : np.nan,
            "B" : np.nan,
            "H" : np.nan,
            "Save" : np.nan
        }

        # sad
        self.alg_dict["sad"][r_id] = {
            "q" : np.nan,
            "n" : np.nan,
            "a0" : np.nan
        }

        # metroman    
        self.alg_dict["metroman"][r_id] = {
             "q" : np.nan,
             "na" : np.nan,
             "x1" : np.nan,
             "a0" : np.nan
        }

        # sic4dvar
        self.alg_dict["sic4dvar"][r_id] = {
            "q31": np.nan,
            "q5": np.nan,
            "n": np.nan,
            "a0": np.nan
        }

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
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if logged:
                return np.exp(np.nanmean(chains, axis=0))
            else:
                return np.nanmean(chains, axis=0)
