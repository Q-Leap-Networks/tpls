#!/bin/bash

# Example batch script for setting up TPLS code on ARCHER

module swap PrgEnv-cray PrgEnv-gnu/5.0.41

module load cray-petsc/3.5.2.1

module load cray-netcdf-hdf5parallel
