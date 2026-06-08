### Efficiency and scalability comparison

This example corresponds to Section 7 in the paper.

Figure 18 (Wall-clock time comparison) is simply obtained by tabulating the wall-clock time for each example. The plotting script is `raw_cost.py`.

For the examples below, similarly to HSGN and SGN, wall-clock time will be needed

For Figure 19 (Wall-clock time vs Number of grid cells), plot the results by running `python mx_alpha_comparison.py` (data reported in the paper already tabulated in the script). To get your own timing comparison:
1. Generate the data for SWE and mbHSGN running `bash run_automated_mx_alpha_refinement.sh`
2. Retrieve the wall-clock times by running `bash mx_alpha_comparison.sh`
One can proceed analogously for SGN inside the directory `SGN`.

For Figures 20 and 21 (Parallel scaling), plot the results by running `python weak_scaling_comparison.py` and `python strong_scaling_comparison.py` (data reported in the paper already tabulated in the scripts). To get your own timing comparison:
1. Generate the data for SWE and mbHSGN running `bash run_automated_weak_scaling_fixed_grid.sh` and `bash run_automated_strong_scaling_fixed_grid.sh`
2. Retrieve the wall-clock times by running `bash weak_scaling_fixed_grid_comparison.sh` and `bash strong_scaling_comparison_fixed_grid.sh`
Similarly for SGN inside the directory `SGN`. One can also generate data with AMR with the scripts `run_automated_weak_scaling.sh`, `run_automated_strong_scaling.sh `and the corresponding comparison scripts.






