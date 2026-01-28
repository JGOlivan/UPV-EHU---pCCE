############################################################################################################
################        This file includes all the Packages needed for the main file        ################
############################################################################################################

using ArgParse, Distributed, MKL, Revise, ProgressMeter, Dates, BenchmarkTools, LinearAlgebra # essential packages
using SlurmClusterManager # only for HPC cluster usage

using Plots
using LsqFit        # only needed for plotting
using Polynomials   # only needed for plotting
using LaTeXStrings  # for writting LaTeX stuff (i think only for plotting)