#!/bin/bash

#Running the experiment in the current directory with different combinations of formal accuracy

#####Defining parameters

#Save current directory
parent_dir=$(pwd)

flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py

export OMP_NUM_THREADS=8
export OMP_STACKSIZE=32M

ulimit -s unlimited

#prefixes
output_prefix="_output_automated"

list_equations=("SWE mbHSGN HSGN")
list_idx_eqs=("0 -1 -2")

#Name of files
makefile="Makefile"
setrun="setrun.py"

#####################################
#####################################
#UPDATE DATA AND RUN SWE SIMULATION
new_line_equations="Bouss_Hyp_rel = 0   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
else
    sed -i "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
fi
make new
#Run simulation
make -B .output "OUTDIR=_output_SWE"

#####################################
#####################################
#UPDATE DATA AND RUN SWE SIMULATION
new_line_equations="Bouss_Hyp_rel = -1   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
else
    sed -i "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
fi
#Run simulation
make -B .output "OUTDIR=_output_mbHSGN"

#####################################
#####################################
#UPDATE DATA AND RUN SWE SIMULATION
new_line_equations="Bouss_Hyp_rel = -2   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
else
    sed -i "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
fi
#Run simulation
make -B .output "OUTDIR=_output_HSGN"
#####################################
#####################################