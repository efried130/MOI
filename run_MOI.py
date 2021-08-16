"""Python script to execute MOI operations.
"""

# Standard imports
from pathlib import Path

# Local imports
from src.Input import Input
from src.Integrate import Integrate
from src.Output import Output

# Constants
INPUT_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Analysis/SWOT/Discharge/Confluence/moi_rundir")
OUTPUT_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Analysis/SWOT/Discharge/Confluence/moi_outdir")

def main():
    input = Input(INPUT_DIR / "flpe", INPUT_DIR / "basin.json", INPUT_DIR / "sos", INPUT_DIR / "swot" )
    input.extract_alg()
    input.extract_swot()

    integrate = Integrate(input.alg_dict, input.basin_dict, input.sos_dict, input.obs_dict)
    integrate.integrate()

    output = Output(input.basin_dict, OUTPUT_DIR, integrate.integ_dict)
    output.write_output()

if __name__ == "__main__":
    from datetime import datetime
    start = datetime.now()
    main()
    end = datetime.now()
    print(f"Execution time: {end - start}")
