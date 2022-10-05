# MOI

MOI serves as the integrator module to the Confluence workflow. It extracts reach-level FLPE algorithm results and SoS data to integrate results on the basin-level. This module writes output to a specified output directory.

# installation

# setup

# execution

**Example Run**
```
%run /Users/mtd/GitHub/SWOT-confluence/moi/run_MOI.py basin.json -v 'unconstrained' 0 
```

So the command line arguments are the basin file, the verbose flag, and the branch name, where branch name can be either constrained or unconstrained. The final argument is the basin number, to be provided only for offline runs.
