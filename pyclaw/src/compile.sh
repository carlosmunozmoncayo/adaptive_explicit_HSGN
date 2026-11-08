#!/usr/bin/bash

rm -rf *.so
#Riemann solvers
python -m numpy.f2py -c linear_solver/lu_solver.f90 rp1_dispersive_SWEs_fwave_splitting.f90 -m rs_bouss_hyprel_fwave_splitting_1D -llapack -lblas
#Total fluctuation
python -m numpy.f2py -c total_fluct.f90 -m dispersive_tfluct
#Source terms
python -m numpy.f2py -c linear_solver/lu_solver.f90 source_terms.f90 -m fortran_source_terms -llapack -lblas
#1D windowed convolution
python -m numpy.f2py -c conv1D.f90 -m conv1D

#Sharpclaw 1D (to use characteristicwise reconstruction)
SHARPCLAW=${CLAW}/pyclaw/src/pyclaw/sharpclaw
python -m numpy.f2py -m bouss_hyprel_sharpclaw1 -c evec.f90 ${SHARPCLAW}/flux1.f90 ${SHARPCLAW}/ClawParams.f90 ${SHARPCLAW}/workspace.f90 ${SHARPCLAW}/weno.f90 ${SHARPCLAW}/reconstruct.f90