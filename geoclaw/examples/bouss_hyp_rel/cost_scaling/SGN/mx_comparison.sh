#!/bin/bash

BASE_DIR="_output_mx_refinement"

# Define the ordered list of proc directories
PROC_ORDER=(proc_4 proc_8 proc_16 proc_32 proc_64 proc_96)
list_mx=("25" "50" "100" "200" "400")
list_alphas=("3.0" "4.2426" "6.0" "8.4853")

echo "=== SGN ==="
for mx in "${list_mx[@]}"; do #Print mx without endline
    echo -n "  mx_$mx "
done
echo "" #Print endline
for mx in "${list_mx[@]}"; do
    timing_file=$BASE_DIR/"mx_$mx/timing.txt"
    if [[ -f "$timing_file" ]]; then
        # Extract the first number from the "Total time:" line
        time_val=$(grep "Total time:" "$timing_file" | awk '{print $3}')
        # echo "  mx_$mx: $time_val"
        echo -n "  $time_val "
    else
        echo -n "  [not found] "
    fi
done
echo "" #Print endline

