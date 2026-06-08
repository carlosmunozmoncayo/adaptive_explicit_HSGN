#!/bin/bash

#Running a 1d mbHSGN (EDC) simulation with a forcing source term

#Save current directory
parent_dir=$(pwd)

#Flags to look for 
flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py

#Name of files
makefile="Makefile"
setrun="setrun.py"
setplot="setplot.py"

#Warning
donot="PLEASE DO NOT DELETE THIS COMMENT"


#####################################
#####################################
#UPDATE DATA AND RUN EDC SIMULATION 1D
#Cd, update setrun.py, run simulation, come back
cd 1d_radial
new_line_equations="Bouss_Hyp_rel = -1   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
if [[ "$OSTYPE" == "darwin"* ]]; then
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
else
    sed -i "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
fi
#Run simulation
make -B .output "OUTDIR=_output_EDC_paper"
cd ${parent_dir}