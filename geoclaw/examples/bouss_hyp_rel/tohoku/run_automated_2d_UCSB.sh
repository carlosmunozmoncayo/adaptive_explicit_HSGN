#!/bin/bash

#####Defining parameters

#Save current directory
parent_dir=$(pwd)

flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py

flag_fujii="flag_fujii"
line_fujii="#    dtopo_data.dtopofiles.append([1,1,1,dtopodir+'Fujii.txydz']) #flag_fujii"
flag_ucsb="flag_ucsb"
line_ucsb="    dtopo_data.dtopofiles.append([1,1,4,dtopodir+'UCSB3.txydz']) #flag_ucsb"

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

#Comment Fujii source and inject UCSB source
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^.*${flag_fujii}.*$/${line_fujii}/" "${setrun}"
    sed -i '' "s/^.*${flag_ucsb}.*$/${line_ucsb}/" "${setrun}"
else
    sed -i "s/^.*${flag_fujii}.*$/${line_fujii}/" "${setrun}"
    sed -i "s/^.*${flag_ucsb}.*$/${line_ucsb}/" "${setrun}"
fi

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
make -B .output "OUTDIR=_output_UCSB_SWE"

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
make -B .output "OUTDIR=_output_UCSB_mbHSGN"

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
make -B .output "OUTDIR=_output_UCSB_HSGN"
#####################################
#####################################