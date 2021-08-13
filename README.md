# prediagnostics

MOI serves as the integrator module to the Confluence workflow. It extracts reach-level FLPE algorithm results and SoS data to integrate results on the basin-level. This module writes output to a specified output directory.

TO DO:
- Implement Input.extract_sos()
- Modify Input.extract_alg() to inlcude MetroMan results, locate HiVDI n parameter, locate MOMMA A0 parameter.
- Integrate.moi_params initialization
- Integrate.stage1_estimate initialization
- Integrate.stage2_estimate addition ??
- Integrate.get_pre_mean_q implementation
- Integrate.integrate implementation
- Output.write_output: Add optional attribute metadata, decide on how to store output of FLPE integration


# installation

# setup

# execution