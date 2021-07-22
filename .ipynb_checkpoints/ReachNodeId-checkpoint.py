# Standard imports
from datetime import datetime
from os import scandir
from pathlib import Path

# Third party imports
from netCDF4 import Dataset

class ReachNodeId:
    """ Extracts reach and node identifiers from SWORD and stores in the SoS.

    Attributes
    ----------
    FILL_VALUE: integer
        Fill value for NetCDF4 missing data.
    sword_dir: Path
        Path to SWORD directory
    sos_dir: Path
        Path to SoS directory
    
    Methods
    -------
    extract_data()
        Extracts SWORD data and stores it in the SoS.
    get_reaches(sword, sos)
        Define dimensions and reach_id variable and write identifier data to sos.
    get_nodes(sword, sos)
        Define dimensions and id variables and write identifier data to sos.
    """

    FILL_VALUE = -9999

    def __init__(self, sword_dir, sos_dir):
        self.sword_dir = sword_dir
        self.sos_dir = sos_dir

    def extract_data(self):
        """Extracts SWORD identifier data and stores it in the SoS matching
        SWORD structure.
        """

        with scandir(self.sword_dir) as entries:
            for sword_file in entries:
                sword = Dataset(Path(sword_file))
                name = sword_file.name.split(".nc")[0] + "_SOS.nc"
                sos = Dataset(self.sos_dir.joinpath(name), 'w', format="NETCDF4")
                sos.Name = sword.Name
                sos.production_date = datetime.now().strftime('%d-%b-%Y %H:%M:%S')
                self.write_reaches(sword, sos)
                self.write_nodes(sword, sos)
                sword.close()
                sos.close()
    
    def write_reaches(self, sword, sos):
        """Write reach_id variable and associated dimension to the SoS."""
        
        sos_reach = sos.createGroup("reaches")
        sos_reach.createDimension("num_reaches", None)
        reach_var = sos_reach.createVariable("reach_id", "i8", ("num_reaches",),
            fill_value=self.FILL_VALUE)
        reach_var.format = sword["reaches"]["reach_id"].format
        reach_var[:] = sword["reaches"]["reach_id"][:]

    def write_nodes(self, sword, sos):
        """Write node_id and reach_id variables with associated dimension to the
        SoS."""

        sos_node = sos.createGroup("nodes")
        sos_node.createDimension("num_nodes", None)
        node_var = sos_node.createVariable("node_id", "i8", ("num_nodes",), 
            fill_value=self.FILL_VALUE)
        node_var.format = sword["nodes"]["node_id"].format
        node_var[:] = sword["nodes"]["node_id"][:]
        reach_var = sos_node.createVariable("reach_id", "i8", ("num_nodes",), 
            fill_value=self.FILL_VALUE)
        reach_var.format = sword["nodes"]["reach_id"].format
        reach_var[:] = sword["nodes"]["reach_id"][:]
    
if __name__ == "__main__":
    #SWORD_DIR = Path("/home/nikki/Documents/confluence/data/sword/SWORD_v07/Reaches_Nodes/netcdf")
    SWORD_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Data/SWOT/SWORD/Reaches_Nodes/netcdf")
    #SOS_DIR = Path("/home/nikki/Documents/confluence/data/sos/sos/netcdf")
    SOS_DIR = Path("/Users/mtd/OneDrive - The Ohio State University/Data/SWOT/SoS/euro_lisflood/SoS/tmp")

    rnid = ReachNodeId(SWORD_DIR, SOS_DIR)
    rnid.extract_data()
