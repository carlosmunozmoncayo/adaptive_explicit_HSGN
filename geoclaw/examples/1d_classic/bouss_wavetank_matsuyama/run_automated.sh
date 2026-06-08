#!/bin/bash

#Save current directory
parent_dir=$(pwd)

flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py
flag_EDC_alpha="flag_EDC_alpha" #Goes in setrun.py
flag_grid_mx="flag_grid_mx" #Goes in make_celledges.py

#Name of files
makefile="Makefile"
setrun="setrun.py"

#Warning
donot="PLEASE DO NOT DELETE THIS COMMENT"


new_line_grid_mx="mx = 8000    #${flag_grid_mx} ${donot}"
#Insert line in make_celledges.py
sed -i '' "s/^.*${flag_grid_mx}.*$/${new_line_grid_mx}/" "${parent_dir}/make_celledges.py"
#Make topography
make -B topo


##############################
#UPDATE DATA AND RUN HSGN SIMULATION
new_line_equations="    rundata.bouss_data.bouss_equations = -2    #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
#Insert new alpha line in setrun.py
new_line_alpha="    alpha = "3.0"           #${flag_EDC_alpha} ${donot}"    
sed -i '' "s/^.*${flag_EDC_alpha}.*$/${new_line_alpha}/" "${setrun}"
#Run simulation
make -B .output OUTDIR="_output_HSGN"
#####################################
#UPDATE DATA AND RUN SWE SIMULATION
new_line_equations="    rundata.bouss_data.bouss_equations = 0   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
#Run simulation
make -B .output OUTDIR="_output_SWE"
#####################################
#UPDATE DATA AND RUN SGN SIMULATION
new_line_equations="    rundata.bouss_data.bouss_equations = 2   #${flag_boussdata_equations} ${donot}"
#Insert line in setrun.py
sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
#Run simulation
make -B .output OUTDIR="_output_SGN"
#####################################
