"""Python script to execute MOI operations.
"""

# Standard imports
import json
import os
from pathlib import Path
import sys

# Local imports
from moi.Input import Input
from moi.Integrate import Integrate
from moi.Output import Output


def get_basin_data(basin_json,index_to_run):
    """Extract reach identifiers and return dictionary.
    
    Dictionary is organized with a key of reach identifier and a value of
    SoS file as a Path object.
    """
    #index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    #index = 0

    if index_to_run == -235:
        index=int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    else:
        index=index_to_run

    with open(basin_json) as json_file:
        data = json.load(json_file)


    return {
        #"basin_id" : int(data[index]["basin_id"]),
        "basin_id" : data[index]["basin_id"], #hope it's ok not to have basin ids always integers?
        "reach_ids" : data[index]["reach_id"],
        "sos" : data[index]["sos"],
        "sword": data[index]["sword"]
    }

def main():

    # verbose 
    try: 
        VerboseFlag=sys.argv[2]
        if VerboseFlag == '-v': Verbose=True
    except IndexError:
        Verbose=False

    #branch
    try:
        Branch=sys.argv[3]
    except IndexError:
        Branch='unconstrained'

    #context
    try:
        index_to_run=int(sys.argv[4]) #integer
    except IndexError:
        index_to_run=-235

    #data directories
    if index_to_run == -235:
        INPUT_DIR = Path("/mnt/data/input")
        FLPE_DIR = Path("/mnt/data/flpe")
        OUTPUT_DIR = Path("/mnt/data/output")
    else:
        INPUT_DIR = Path("/Users/mtd/Analysis/SWOT/Discharge/Confluence/paper_debug/mnt/input")
        FLPE_DIR = Path("/Users/mtd/Analysis/SWOT/Discharge/Confluence/paper_debug/mnt/flpe")
        OUTPUT_DIR = Path("/Users/mtd/Analysis/SWOT/Discharge/Confluence/paper_debug/moi_outputs_dev")

    #basin data
    try:
        basin_json = INPUT_DIR.joinpath(sys.argv[1])
    except IndexError:
        basin_json = INPUT_DIR.joinpath("basin.json") 


    basin_data = get_basin_data(basin_json,index_to_run)

    print('Running ',Branch,' branch.')

    input = Input(FLPE_DIR, INPUT_DIR / "sos/", INPUT_DIR / "swot", INPUT_DIR / "sword", basin_data)

    input.extract_sos()
    input.extract_alg()
    input.extract_swot()
    input.extract_sword()

    integrate = Integrate(input.alg_dict, input.basin_dict, input.sos_dict, input.sword_dict,input.obs_dict,Branch,Verbose)
    integrate.integrate()

    output = Output(input.basin_dict, OUTPUT_DIR, integrate.integ_dict, integrate.alg_dict, integrate.obs_dict, input.sword_dir)
    output.write_output()
    output.write_sword_output(Branch)

if __name__ == "__main__":
    from datetime import datetime
    start = datetime.now()
    main()
    end = datetime.now()
    print(f"Execution time: {end - start}")
