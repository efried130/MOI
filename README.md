# prediagnostics

MOI serves as the integrator module to the Confluence workflow. It extracts reach-level FLPE algorithm results and SoS data to integrate results on the basin-level. This module writes output to a specified output directory.

TO DO:
- Basin JSON file that identifies basin and corresponding reach identifiers and SoS file
- Implement Input.extract_sos()
- Modify Input.extract_alg() to extract basin-level data instead of reach-level data, locate HiVDI n parameter, locate MOMMA A0 parameter.
- Input.__get_ids(basin_json): organize with key of basin number and value of dict with reach ids and sos file.
- Integrate.moi_params initialization
- Integrate.stage1_estimate initialization
- Integrate.stage2_estimate addition ??
- Integrate.integrate implementation
- Output.write_output: Add optional attribute metadata, decide on how to store output of FLPE integration


# installation

# setup

# execution