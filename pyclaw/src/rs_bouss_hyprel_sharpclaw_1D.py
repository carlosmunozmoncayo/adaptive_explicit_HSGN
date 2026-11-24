"""Riemann solver for PyClaw's SharpClawSolver1D (EDC/mbHSGN, BBBD/HSGN, DR):

    import rs_bouss_hyprel_sharpclaw_1D as rs_sharpclaw
    solver = pyclaw.SharpClawSolver1D(rs_sharpclaw)

It is the compiled module rs_bouss_hyprel_fwave_splitting_1D (same cparam common block), with the
entry point rp1_sharpclaw used as rp1: SharpClaw reconstructs h (b is piecewise constant), so the
non-conservative products of EDC/BBBD use the jump of h instead of the jump of eta. Use it together
with sharpclaw_source_step_EDC / sharpclaw_source_step_BBBD from source_terms.py.
The classic solver uses rs_bouss_hyprel_fwave_splitting_1D directly.
"""
import rs_bouss_hyprel_fwave_splitting_1D as _rs

rp1 = _rs.rp1_sharpclaw
cparam = _rs.cparam
