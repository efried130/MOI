"""Python script to execute MOI operations.
"""

# Standard imports
import json
import os
from pathlib import Path
import sys

# Local imports
from src.Input import Input
from src.Integrate import Integrate
from src.Output import Output

# Constants
INPUT_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Analysis/SWOT/Discharge/Confluence/moi_rundir")
FLPE_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Analysis/SWOT/Discharge/Confluence/moi_rundir/flpe")
OUTPUT_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Analysis/SWOT/Discharge/Confluence/moi_outdir")

def get_basin_data(basin_json):
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

def main():

    try:
        basin_json = INPUT_DIR.joinpath(sys.argv[1])
    except IndexError:
        basin_json = INPUT_DIR.joinpath("basin.json") 
    basin_data = get_basin_data(basin_json)

    input = Input(FLPE_DIR, INPUT_DIR / "sos", INPUT_DIR / "swot", basin_data)
    input.extract_alg()
    input.extract_swot()

    integrate = Integrate(input.alg_dict, input.basin_dict, input.sos_dict, input.obs_dict)
    integrate.integrate()

    output = Output(input.basin_dict, OUTPUT_DIR, integrate.integ_dict, integrate.alg_dict, integrate.obs_dict)
    output.write_output()

if __name__ == "__main__":
    from datetime import datetime
    start = datetime.now()
    main()
    end = datetime.now()
    print(f"Execution time: {end - start}")
