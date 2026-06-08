#!/bin/bash

#Running the experiment in the current directory with different combinations of formal accuracy

#####Defining parameters

#Save current directory
parent_dir=$(pwd)

flag_boussdata_equations="flag_bouss_data.bouss_equations" #Goes in setrun.py


#prefixes
output_prefix="_output_mx_refinement"

list_final_times=("150."  "150." "150." "150." "150.")
list_dom_sizes=("5000.0" "5000.0" "5000.0" "5000.0" "5000.0")
# num_processors=("4" "8" "16" "32" "64" "96")

list_mx=("25" "50" "100" "200" "400")
list_my=("25" "50" "100" "200" "400")


flag_dom_sizex="flag_dom_sizex" #Goes in setrun.py
flag_dom_sizey="flag_dom_sizey" #Goes in setrun.py
flag_final_time="flag_final_time" #Goes in setrun.py
flag_mx="flag_mx" #Goes in setrun.py
flag_my="flag_my" #Goes in setrun.py

donotdelete="PLEASE DO NOT DELETE THIS COMMENT" #String to identify lines to be replaced in setrun.py
#Name of files
makefile="Makefile"
setrun="setrun.py"



export OMP_STACKSIZE=32M
if [[ "$OSTYPE" == "darwin"* ]]; then
    ulimit -s hard
else
    ulimit -s unlimited
fi
#Update number of processors
export FFLAGS='-O2 -fopenmp'
export OMP_NUM_THREADS=8
export BOUSS_MPI_PROCS=8
make new

for i in "${!list_final_times[@]}"; do
    tf="${list_final_times[$i]}"
    L="${list_dom_sizes[$i]}"
    np="8"
    mx="${list_mx[$i]}"
    my="${list_my[$i]}"

    echo "Running simulations with tfinal = ${tf}, L = ${L}, np = ${np}, mx = ${mx} and my = ${my}">>log_automated.txt

    #UPDATE DATA AND RUN SWE SIMULATION
    new_line_dom_size0="    clawdata.upper[0] = ${L}       #${flag_dom_sizex}  ${donotdelete}"
    new_line_dom_size1="    clawdata.upper[1] = ${L}       #${flag_dom_sizey}  ${donotdelete}"
    new_line_final_time="        clawdata.tfinal = ${tf}  #${flag_final_time}  ${donotdelete}"
    new_line_mx="    clawdata.num_cells[0] = ${mx}      #${flag_mx}  ${donotdelete}"
    new_line_my="    clawdata.num_cells[1] = ${my}      #${flag_my}  ${donotdelete}"
    
    #Insert line in setrun.py
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/^.*${flag_dom_sizex}.*$/${new_line_dom_size0}/" "${setrun}"
        sed -i '' "s/^.*${flag_dom_sizey}.*$/${new_line_dom_size1}/" "${setrun}"
        sed -i '' "s/^.*${flag_final_time}.*$/${new_line_final_time}/" "${setrun}"
        sed -i '' "s/^.*${flag_mx}.*$/${new_line_mx}/" "${setrun}"
        sed -i '' "s/^.*${flag_my}.*$/${new_line_my}/" "${setrun}"
    else
        sed -i "s/^.*${flag_dom_sizex}.*$/${new_line_dom_size0}/" "${setrun}"
        sed -i "s/^.*${flag_dom_sizey}.*$/${new_line_dom_size1}/" "${setrun}"
        sed -i "s/^.*${flag_final_time}.*$/${new_line_final_time}/" "${setrun}"
        sed -i "s/^.*${flag_mx}.*$/${new_line_mx}/" "${setrun}"
        sed -i "s/^.*${flag_my}.*$/${new_line_my}/" "${setrun}"
    fi
    
    #Run simulation
    make data
    echo "Running SGN simulation with tfinal = ${tf}, L = ${L} and np = ${np}">>log_automated.txt
    make -B .output "OUTDIR=${output_prefix}/mx_${mx}"

done