#!/bin/sh
# Adios implementation: Ebru & Matthieu
# Princeton, August 2013

if [ ! -f ../constants.h ]; then
        echo Error: add constants.h
        exit
fi
if [ ! -f ../values_from_mesher.h ]; then
        echo Error: add values_from_mesher.h 
        exit
fi

if [ ! -f ../precision.h ]; then 
	echo Error: add precision.h 
	exit
fi 

adios_link=`adios_config -lf`
adios_inc=`adios_config -cf`

ftn -O3 -c adios_helpers_definitions.f90 $adios_inc
ftn -O3 -c adios_helpers_writers.f90 $adios_inc
ftn -O3 -c adios_helpers.f90 $adios_inc

ftn -O0 -g -fcheck=all -o xcompute_cg_direction compute_cg_direction_single_beta.f90 exit_mpi.f90 adios_helpers.o adios_helpers_definitions.o adios_helpers_writers.o $adios_inc $adios_link
#mpif90 -O3 -o xcompute_cg_direction compute_cg_direction_single_beta.f90 exit_mpi.f90
