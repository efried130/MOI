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

    def __init__(self, alg_dir, sos_dir, swot_dir, sword_dir,basin_data):
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
        self.sword_dir = sword_dir
        self.swot_dir = swot_dir

    def extract_sos(self):
        """Extracts and stores SoS data in sos_dict.
        
        Parameters
        ----------
        """

        sosfile=self.sos_dir.joinpath(self.sos_dir, self.basin_dict['sos'])

        sos_dataset=Dataset(sosfile)
    
        sosreachids=sos_dataset["reaches/reach_id"][:]
        sosQbars=sos_dataset["model/mean_q"][:]

        self.sos_dict={}
        for reach in self.basin_dict['reach_ids']:
            self.sos_dict[reach]={}
            k=np.argwhere(sosreachids==np.int64(reach))
            k=k[0,0]
            self.sos_dict[reach]['Qbar']=sosQbars[k]

    def extract_sword(self):
        """Extracts and stores SWORD data in sword_dict.
        
        Parameters
        ----------
        """
        swordfile=self.sword_dir.joinpath(self.sword_dir, self.basin_dict['sword'])
        sword_dataset=Dataset(swordfile)

        self.sword_dict={} #organized by field rather than by reaches

        # grab sizes of the data
        dimfields=['orbits','num_domains','num_reaches']
        for field in dimfields:
            self.sword_dict[field]=sword_dataset['reaches'].dimensions[field].size    

        # grab data    
        reachfields=['reach_id','facc','n_rch_up','n_rch_down','rch_id_up','rch_id_dn','swot_obs','swot_orbits']
        for field in reachfields:
            self.sword_dict[field]=sword_dataset['reaches/' + field][:]


    def extract_swot(self):
 
        self.obs_dict={}

        for reach in self.basin_dict['reach_ids']:
             swotfile=self.swot_dir.joinpath(reach+'_SWOT.nc')
             swot_dataset = Dataset(swotfile)

             self.obs_dict[reach]={}
             nt = swot_dataset.dimensions['nt'].size
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

            if not mm_file:
                mm_file=Path('dir/that/does/not/exist')  #this sets mm_file.exists() to false
            else: 
                mm_file = Path(mm_file[0]) 

            self.__extract_valid(r_id, gb_file, hv_file, mo_file, sd_file, mm_file, sv_file)

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
        if gb_file.exists():
            gb = Dataset(gb_file, 'r', format="NETCDF4")
            self.alg_dict["geobam"][r_id] = {
                "s1-flpe-exists": True,
                "q": np.array(self.__get_gb_data(gb, "q", "q", True)),
                "n": np.array(self.__get_gb_data(gb, "logn", "mean", True)),
                "a0": 1.0    # TODO temp value until work out neoBAM A0
            }
            gb.close()
        else:
            self.alg_dict["geobam"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "n" : np.nan,
                "a0" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

        # hivdi
        if hv_file.exists():
            hv = Dataset(hv_file, 'r', format="NETCDF4")
            self.alg_dict["hivdi"][r_id] = {
                "s1-flpe-exists": True,
                "q" : hv["reach"]["Q"][:].filled(np.nan),
                "alpha" : hv["reach"]["alpha"][:].filled(np.nan),  
                "beta" : hv["reach"]["beta"][:].filled(np.nan),  
                "a0" : hv["reach"]["A0"][:].filled(np.nan)
            }
            hv.close()
        else:
            self.alg_dict["hivdi"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "alpha" : np.nan,
                "beta" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

        # momma
        if mo_file.exists():
            mo = Dataset(mo_file, 'r', format="NETCDF4")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                self.alg_dict["momma"][r_id] = {
                    "s1-flpe-exists": True,
                    "q" : mo["Q"][:].filled(np.nan),
                    "B" : mo["zero_flow_stage"][:].filled(np.nan),
                    "H" : mo["bankfull_stage"][:].filled(np.nan),                                  
                    "Save" : np.nanmean(mo["slope"][:].filled(np.nan))
                }
            mo.close()
        else:
            self.alg_dict["momma"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "B" : np.nan,
                "H" : np.nan,
                "Save" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

        # sad
        if sd_file.exists():
            sd = Dataset(sd_file, 'r', format="NETCDF4")
            self.alg_dict["sad"][r_id] = {
                "s1-flpe-exists": True,
                "q" : sd["Qa"][:].filled(np.nan),
                "n" : sd["n"][:].filled(np.nan),
                "a0" : sd["A0"][:].filled(np.nan)
            }
            sd.close()
        else:
            self.alg_dict["sad"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "n" : np.nan,
                "a0" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

        # metroman    
        if mm_file.exists():
            mm = Dataset(mm_file, 'r', format="NETCDF4")
            index = np.where(mm["reach_id"][:] == int(r_id))
            self.alg_dict["metroman"][r_id] = {
                 "s1-flpe-exists": True,
                 "q" : mm["allq"][index].filled(np.nan),
                 "na" : mm["nahat"][index].filled(np.nan),
                 "x1" : mm["x1hat"][index].filled(np.nan),
                 "a0" : mm["A0hat"][index].filled(np.nan)
            }
            mm.close()
        else:
            self.alg_dict["metroman"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "na" : np.nan,
                "x1" : np.nan,
                "a0" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

        # sic4dvar
        if sv_file.exists():
            sv = Dataset(sv_file, 'r', format="NETCDF4")
            self.alg_dict["sic4dvar"][r_id] = {
                #"q31": sv["Qalgo31"][:].filled(np.nan),#unclear which of these to use
                "s1-flpe-exists": True,
                "q": sv["Qalgo31"][:].filled(np.nan),
                "q5": sv["Qalgo5"][:].filled(np.nan),
                "n": sv["n"][:].filled(np.nan),
                "a0": sv["A0"][:].filled(np.nan)
            }
            sv.close()
        else:
            self.alg_dict["sic4dvar"][r_id] = { 
                "s1-flpe-exists" : False ,
                "q" : np.nan,
                "q5" : np.nan,
                "n" : np.nan,
                "a0" : np.nan,
                "qbar" : self.sos_dict[r_id]['Qbar']
            }

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

    def __get_gb_data(self, gb, group, pre, logged):
        """Return geoBAM data as a numpy array.
        
        Parameters
        ----------
        gb: netCDF4.Dataset
            NetCDF file dataset to extract discharge time series
        group: str
            string name of group to access chains
        pre: str
            string prefix of variable name
        logged: bool
            boolean indicating if result is logged
        """

        chain1 = gb[group][f"{pre}1"][:].filled(np.nan)
        chain2 = gb[group][f"{pre}2"][:].filled(np.nan)
        chain3 = gb[group][f"{pre}3"][:].filled(np.nan)
        chains = np.vstack((chain1, chain2, chain3))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if logged:
                return np.exp(np.nanmean(chains, axis=0))
            else:
                return np.nanmean(chains, axis=0)
