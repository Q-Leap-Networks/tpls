#!/bin/bash

# Example batch script for setting up TPLS code on ARCHER

module swap PrgEnv-cray PrgEnv-gnu

module load cray-petsc/3.6.1.0

module load cray-netcdf-hdf5parallel
