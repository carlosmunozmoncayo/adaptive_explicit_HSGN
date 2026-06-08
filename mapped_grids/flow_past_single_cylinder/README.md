### Immersed cylindrical obstacle

This example corresponds to Section 6.3 of the paper.
While it is not based in GeoClaw, it uses a signficant portion of the underlying AMRClaw code.
To generate the data for mbHSGN, run the following commands:
1. make data
2. make .output
For the SWE, take `alpha_dimensionless_relaxation=0`, turn off transverse propagation: `clawdata.transverse_waves = 'none'`, and reduce the CFL number:`clawdata.cfl_desired = 0.45` in `setrun.py`. To compensate for the decreased CFL number, also increase the initial refinement ratio to 4 (instead of 2) in `setrun.py`.

To generate Figures 14 and 15 in the paper, run:
1. python plots_surface_paper.py (make sure the directory Figures exists)
