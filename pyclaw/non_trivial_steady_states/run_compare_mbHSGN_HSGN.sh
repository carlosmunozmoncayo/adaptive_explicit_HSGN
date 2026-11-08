#!/bin/bash

python simulation_comparison_steady_SGN.py outdir=_output_HSGN
python simulation_comparison_steady_SGN.py outdir=_output_mbHSGN bouss_equations=EDC
python compare_mbHSN_HSGN.py
