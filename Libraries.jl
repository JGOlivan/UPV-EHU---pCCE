############################################################################################################
################        This file includes all the Packages needed for the main file        ################
############################################################################################################

using ArgParse, Distributed, MKL, Revise, ProgressMeter, Dates, BenchmarkTools, LinearAlgebra # essential packages
using SlurmClusterManager # only for HPC cluster usage

using Plots
using LaTeXStrings