#!/bin/bash

#Set some environment variables for the HSGN_MAPPED project
#First check if $CLAW is set
if [ -z "$CLAW" ]; then
    echo "Error: \$CLAW is not set. Please set it to the root directory of your clawpack installation."
    exit 1
fi

#Sets some environment variables, run with source
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HSGN_MAPPED_SRC="$(cd ${script_dir}/src && pwd)"
export HSGN_MAPPED_SRC="$HSGN_MAPPED_SRC"
echo "HSGN_MAPPED_SRC is set to $HSGN_MAPPED_SRC"