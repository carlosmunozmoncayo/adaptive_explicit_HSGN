### Tohoku tsunami

Setup taken from https://github.com/rjleveque/ImplicitAMR-paper. Please refer to that repository for instructions on how to generate the data (directory `tohoku_sgn_hawaii`) for the SGN system (requires PETSc), and for details on the sources and gauge data.

To generate the Figures of Section 6.4 in the paper, run the following commands from the current directory:
1. bash run_automated_2d_Fujii.sh
2. bash run_automated_2d_UCSB.sh
3. python plots_surface_paper.py (make sure the directory Figures exists)
4. python compare_gauges_dart.py

Before plotting, data for SGN (`_output_UCSB_SGN` and `_output_fujii_SGN`) has to be generated.
