############################################################################################################
################                This file installs all the necessary Packages               ################
############################################################################################################

Version = true

import Pkg

function installPkg(Pkgs :: Vector{String})

    """ This function install all the necessary packages for the simulation. If already installed it will do nothing. """
                
    for package in Pkgs if isnothing(Base.find_package(package)) Pkg.add(package); end; end

end

function installPkg(Pkgs :: Dict{String, String})

    """ This function install all the necessary packages for the simulation. If already installed it will do nothing. """
                
    for (package, V) in Pkgs if isnothing(Base.find_package(package)) Pkg.add(Pkg.PackageSpec(;name = package, version = V)); end; end

end

if !Version

    # For installing WITHOUT a particular version

    Pkgs = [
                    "ArgParse",                 
                    "Distributed",              
                    "MKL",                      
                    "Revise",
                    "ProgressMeter",
                    "Dates",
                    "BenchmarkTools",
                    "SlurmClusterManager",

                    "StaticArrays",
                    "StructArrays",
                    "SparseArrays",
                    "LinearAlgebra",
                    "Statistics",
                    "Random",
                    "Combinatorics",
                    "FastExpm",
                    "JuMP",
                    "GLPK"
                ]

else

    # For installation WITH a particular version

    Pkgs = Dict(
                    "ArgParse" => "1.2.0",
                    "Distributed" => "1.11.0",
                    "MKL" => "0.9.0",
                    "Revise" => "3.13.1",
                    "ProgressMeter" => "1.11.0",
                    "Dates" => "1.11.0",
                    "BenchmarkTools" => "1.6.3",
                    "SlurmClusterManager" => "1.1.0",

                    "StaticArrays" => "1.9.16",
                    "StructArrays" => "0.7.2",
                    "SparseArrays" => "1.12.0",
                    "LinearAlgebra" => "1.12.0",
                    "Statistics" => "1.11.1",
                    "Random" => "1.11.0",
                    "Combinatorics" => "1.1.0",
                    "FastExpm" => "1.1.0",
                    "JuMP" => "1.29.3",
                    "GLPK" => "1.2.1"
    )

end

installPkg(Pkgs)