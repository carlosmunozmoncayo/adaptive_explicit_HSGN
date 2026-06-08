#!/bin/bash

#Running the experiment in the current directory with different resulting Froude numbers.

#Save current directory
parent_dir=$(pwd)

flag_outputdir="flag_output_directory" #Goes in Makefile
flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py
flag_EDC_alpha="flag_EDC_alpha" #Goes in setrun.py
flag_froude_number="flag_froude_number" #Goes in make_celledges.py

#prefixes
output_prefix="_output_automated"


list_equations=("SWE EDC SGN")
list_froude_number=("1.104d0" "1.1702d0")
list_EDC_alpha=("3.0")

#Name of files
makefile="Makefile"
setrun="setrun.py"

#Warning
donot="PLEASE DO NOT DELETE THIS COMMENT"


for froude_number in ${list_froude_number[@]}
do
    new_line_froude="    Froude = ${froude_number} !${flag_froude_number} ${donot}"
    #Insert line in qinit.f90
    sed -i '' "s/^.*${flag_froude_number}.*$/${new_line_froude}/" "qinit.f90"

    #Make topography
    make -B topo

    ##############################
    #UPDATE DATA AND RUN EDC SIMULATION
    new_line_equations="    rundata.bouss_data.bouss_equations = -1    # 0=SWE, 1=MS, 2=SGN, -1=EDC    #${flag_boussdata_equations} ${donot}"
    #Insert line in setrun.py
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
    for EDC_alpha in ${list_EDC_alpha[@]}
    do
        new_line_dir="OUTDIR = "${output_prefix}/Froude_${froude_number}/EDC_alpha_${EDC_alpha}"     #${flag_outputdir} ${donot}"
        #Insert line in Makefile
        sed -i '' "s|^.*${flag_outputdir}.*$|${new_line_dir}|" "${makefile}"
        #Insert new alpha line in setrun.py
        new_line_alpha="    alpha = "${EDC_alpha}"           #${flag_EDC_alpha} ${donot}"    
        sed -i '' "s/^.*${flag_EDC_alpha}.*$/${new_line_alpha}/" "${setrun}"
        #Run simulation
        make .output
    done
    #####################################
    #UPDATE DATA AND RUN SGN SIMULATION
    new_line_equations="    rundata.bouss_data.bouss_equations = 2    # 0=SWE, 1=MS, 2=SGN, -1=EDC    #${flag_boussdata_equations} ${donot}"
    #Insert line in setrun.py
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
    new_line_dir="OUTDIR = "${output_prefix}/Froude_${froude_number}/SGN"     #${flag_outputdir} ${donot}"
    #Insert line in Makefile
    sed -i '' "s|^.*${flag_outputdir}.*$|${new_line_dir}|" "${makefile}"
    #Run simulation
    make .output
    #####################################
    #UPDATE DATA AND RUN SWE SIMULATION
    new_line_equations="    rundata.bouss_data.bouss_equations = 0    # 0=SWE, 1=MS, 2=SGN, -1=EDC    #${flag_boussdata_equations} ${donot}"
    #Insert line in setrun.py
    sed -i '' "s/^.*${flag_boussdata_equations}.*$/${new_line_equations}/" "${setrun}"
    new_line_dir="OUTDIR = "${output_prefix}/Froude_${froude_number}/SWE"     #${flag_outputdir} ${donot}"
    #Insert line in Makefile
    sed -i '' "s|^.*${flag_outputdir}.*$|${new_line_dir}|" "${makefile}"
    #Run simulation
    make .output
    #####################################
done
