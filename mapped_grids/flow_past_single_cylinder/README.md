### Immersed cylindrical obstacle

This example corresponds to Section 6.3 of the paper.
While it is not based in GeoClaw, it uses a signficant portion of the underlying AMRClaw code.
To generate the data for mbHSGN, run the following commands:
1. make new
2. make data
3. make .output

For the SWE, take `alpha_dimensionless_relaxation=0` in `setrun.py` and run the commands above with `make .output OUTDIR=_output_SWE`.

To generate Figures 14 and 15 in the paper, run:
1. make mapping
2. python plot_coarse_mapped_grid.py (make sure the directory Figures exists)
3. python plotting_results_pyvista.py
4. python plot_gauge_analysis.py

Notes (reproducibility check, Nov 2026): create the `Figures` directory before plotting; on a headless
machine run the pyvista script with `PYVISTA_OFF_SCREEN=true python plotting_results_pyvista.py`.