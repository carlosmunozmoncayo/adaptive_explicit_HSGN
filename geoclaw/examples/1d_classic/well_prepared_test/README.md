###  Gaussian hump with/without well-prepared initial conditions
This example correponds to Section 5.1 of the paper.

Run `python generate_well_prepared_IC.py` if the output of run with three steps with SGN is available in the directory `_output_WP_data`. Otherwise, well-prepared initial data is already available in `w0_p0.txt`.

To generate the data, set 
`probdata.add_param('wellprepared', 0,  'Well-prepared initial conditions? 0=no, 1=yes')`
in `setrun.py` and run the following commands:
1. `make topo`
2. `make data`
3. `make .output OUTDIR=_output_mbHSGN_alpha3_noWP`

Then, set `probdata.add_param('wellprepared', 1,  'Well-prepared initial conditions? 0=no, 1=yes')` and run the following commands:
1. `make topo`
2. `make data`
3. `make .output OUTDIR=_output_mbHSGN_alpha3_WP`

To generate Figures 5 and 6, run: `python plot_eta_pressure.py`.
To generate Figure 7, run: `python plot_Fourier.py`.


The plotting scripts also need an SGN reference solution in `_output_SGN`:
set `rundata.bouss_data.bouss_equations = 2` in `setrun.py`, then `make data` and `make .output OUTDIR=_output_SGN`
(the 1D Boussinesq library in `geoclaw/src/1d_classic/bouss` uses the SGN parameter alpha=1, as stated in the paper).
Create the `Figures` directory before plotting.

`_output_WP_data` (input of `generate_well_prepared_IC.py`) is obtained with the same SGN setting together with
`clawdata.output_style = 3` (3 steps, one output per step) and `clawdata.dt_variable = False` (dt = dt_initial = 1e-5).
Note that `generate_well_prepared_IC.py` overwrites `w0_p0.txt`.
