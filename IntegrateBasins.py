#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 14:26:17 2021

@author: mtd
"""

from os import scandir
from pathlib import Path
import numpy as np
from datetime import datetime

# Third party import
from netCDF4 import Dataset  

class IntegrateBasins:
    
    FILL_VALUE = -9999

    def __init__(self,sword_dir,out_dir,MOIparams,Stage1Estimates):
        self.sword_dir=sword_dir
        self.out_dir=out_dir
        self.MOIparams=MOIparams
        self.Stage1Estimates=Stage1Estimates
        
    def Integrate(self):

        with scandir(self.sword_dir) as entries:
            for sword_file in entries:
                
                # get sword data for this file 
                sword = Dataset(Path(sword_file))  
                name=sword.Name
                reach_var = sword["reaches"]["reach_id"][:]      
                
                # run integrator on L2 basins
                # the output is new sets of mean flow and flow law parameters
                nreach=np.size(reach_var)
                nvars=2
                Stage2Estimates={
                        "FLP":np.empty([nvars,len(reach_var)]),
                        "Qmean":np.empty([np.size(reach_var)]) }                
                
                Stage2Estimates["FLP"].fill(self.FILL_VALUE)
                Stage2Estimates["Qmean"].fill(self.FILL_VALUE)                
                
                #write result to new SOS netcdf file
                name = sword_file.name.split(".nc")[0] + "_OUT.nc"
                out = Dataset(self.out_dir.joinpath(name), 'w', format="NETCDF4")
                out.Name = sword.Name
                out.production_date = datetime.now().strftime('%d-%b-%Y %H:%M:%S')
                self.write_reaches(sword, out,Stage2Estimates)
                
                sword.close()
                out.close()

        return 
    
    
    def write_reaches(self, sword, out,Stage2Estimates):
        """Write reach_id variable and associated dimension to the SoS."""
        
        out_reach = out.createGroup("reaches")
        out_reach.createDimension("num_reaches", None)
        out_reach.createDimension("num_flps", None)
        
        #reach id
        reach_var = out_reach.createVariable("reach_id", "i8", ("num_reaches",),
            fill_value=self.FILL_VALUE)
        reach_var.format = sword["reaches"]["reach_id"].format
        reach_var[:] = sword["reaches"]["reach_id"][:]
        
        #Qmean
        reach_Qvar = out_reach.createVariable("Qmean", "f4", ("num_reaches",),
            fill_value=self.FILL_VALUE)
        reach_Qvar[:] = Stage2Estimates["Qmean"][:]
        
        #FLPs
        reach_FLPvar = out_reach.createVariable("FLP", "f4", ("num_flps","num_reaches"),
            fill_value=self.FILL_VALUE)
        reach_FLPvar[:] = Stage2Estimates["FLP"][:]        

        
        
        


