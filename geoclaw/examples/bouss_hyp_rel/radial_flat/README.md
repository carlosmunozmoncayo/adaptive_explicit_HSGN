### Radial hump over flat bathymetry

The starting point for this setup was the example in the GeoClaw repository: https://github.com/clawpack/geoclaw/tree/master/examples/bouss/radial_flat

To generate the Figures of Section 5.5 in the paper, run the following commands from the current directory:
1. bash run_automated_1d.sh
2. bash run_automated_2d.sh
3. python plots_surface_paper.py (make sure the directory Figures exists)

To generate the side by side comparison of Figure 1:
1. Repeat the above steps with a fixed grid size of 400x400.
2. Generate the data for SGN (_output_SGN) running the radial_flat example in the GeoClaw repository with the same grid size.
3. Generate the data for Linear Free Surface Euler with the spectral method in linear_water_waves_problem.ipynb.
4. Finally, the comparison can be made with the notebook compare_all_w_Euler.ipynb.

