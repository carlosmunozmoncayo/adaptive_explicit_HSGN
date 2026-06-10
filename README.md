#### Reproducibility repository for the paper: "Adaptive, Efficient, and Scalable Water Wave Modeling with Dispersive Hyperbolic Systems" by Carlos Muñoz-Moncayo and David I. Ketcheson

Inside each directory, there are instructions to set up the required environment variables and reproduce the results in the manuscript. In particular, code corresponding to the different sections can be found in:
- [Section 5.1](geoclaw/examples/1d_classic/well_prepared_test)
- [Section 5.2](pyclaw/solitary_flat)
- [Section 5.3](pyclaw/non_trivial_steady_states)
- [Section 5.4](geoclaw/examples/bouss_hyp_rel/diagonal_head_on_collision)
- [Section 5.5](geoclaw/examples/bouss_hyp_rel/radial_flat)
- [Section 6.1](geoclaw/examples/1d_classic/dam_break_undular_bores)
- [Section 6.2](geoclaw/examples/1d_classic/bouss_wavetank_matsuyama)
- [Section 6.3](mapped_grids/flow_past_single_cylinder)
- [Section 6.4](geoclaw/examples/bouss_hyp_rel/tohoku)
- [Section 7](geoclaw/examples/bouss_hyp_rel/cost_scaling)

[Clawpack](https://www.clawpack.org/) is required to run the code (tested with Clawpack 5.10.0 and 5.13.1). Please make sure that the environment variable CLAW points to your Clawpack installation.

The code has been run on a Linux (Ubuntu) machine with an AMD Ryzen Threadripper 7980X CPU, and a MacBook Pro M2.


