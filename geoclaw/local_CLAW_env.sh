#!/bin/bash

#Sets some environment variables, run with source
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#Local paths
export LOCALGEOCLAW="$(cd ${script_dir} && pwd)"
export LOCALBOUSSLIB="$(cd ${script_dir}/src/1d_classic/bouss && pwd)"
export LOCALBOUSSLIB2D="$(cd ${script_dir}/src/2d/hyp_bouss && pwd)"

#Paths pointing to clawpack (environment variable CLAW should be set to the directory where clawpack is installed)
export LOCALGEOLIB="${CLAW}/geoclaw/src/1d_classic/shallow"
export LOCALGEOLIB2D="${CLAW}/geoclaw/src/2d/shallow"