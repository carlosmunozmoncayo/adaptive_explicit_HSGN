#!/bin/bash

BASE_DIR="_output_weak_scaling_fixed_grid"

# Define the ordered list of proc directories
PROC_ORDER=(proc_1 proc_2 proc_4 proc_8 proc_16 proc_32 proc_64 proc_80 proc_96)

# Loop over each system directory
# for system_dir in "$BASE_DIR"/*/; do
    # system=$(basename "$system_dir")
    # echo "=== $system ==="

    for proc in "${PROC_ORDER[@]}"; do
        timing_file="$BASE_DIR/$proc/timing.txt"

        if [[ -f "$timing_file" ]]; then
            # Extract the first number from the "Total time:" line
            time_val=$(grep "Total time:" "$timing_file" | awk '{print $3}')
            echo "  $proc: $time_val"
        else
            echo "  $proc: [not found]"
        fi
    # done

    echo ""
done