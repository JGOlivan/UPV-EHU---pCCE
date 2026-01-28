#############################################################################################################
################  This file contains all the functions needed for setting up the simulation  ################
#############################################################################################################

function parse_command_line()

    """
        This function parses the input arguments from the command line. Input arguments from the command line (but not the default values) are priority and will overwrite the Input.jl file. This means: if input arguments are given in the command line, the input file is useless.
    """

  s = ArgParseSettings()

  @add_arg_table s begin

    "--Species"
    arg_type    = String
    help        = "The chemical species of the spin bath. It is a String object. It can be a nuclear spin (N14, N15, C13, ...) or an electron spin (e). Defult is N15"
    default     = "N15"
    required    = false

    "--NumberSpins"
    arg_type    = Int64
    help        = "The number of bath spins included in the full quantum mechanical simulation. It is an Integer. Default is 180."
    default     = 180
    required    = false

    "--concentration"
    arg_type    = Float64
    help        = "The concentration of the bath spins in the lattice in parts per million (ppm). It is a Real. Default is 5ppm."
    default     = 5
    required    = false

    "--nPulses"
    arg_type    = Int64
    help        = "The number of pi pulses in the Dynamical Decoupling sequence applied to the central spin. Default is 1. If set to zero, FID is studied."
    default     = 1
    required    = false
    
    "--Protocol"
    arg_type    = String
    help        = "The protocol considered. It is a String. Options are: Free; OneToneRes, TwoToneRes, FullToneRes; OneToneLG, TwoToneLG, FullToneLG; TwoToneHyb. Default is Free."
    default     = "Free"
    required    = false

    "--Omega"
    arg_type    = Float64
    help        = "The Rabi frequency of the bath driving (in MHZ, per tone). Default is 0 (i.e., no bath driving applied)."
    default     = 0
    required    = false

    "--tauMax"
    arg_type    = Float64
    help        = "The duration of a single free evolution period (in μs). Default is 100."
    default     = 100
    required    = false

    "--numSpatialAv"
    arg_type    = Int64
    help        = "The number of spatial configurations considered for averaging, i.e., the number of central spins in the ensemble. Set to 1 for single central spin. Default is 50."
    default     = 50
    required    = false

    "--meanField"
    arg_type    = Bool
    help        = "Set to true to include mean-field averaging in the simulation. It is a Bool. Default is true."
    default     = true
    required    = false

    "--numExtAv"
    arg_type    = Int64
    help        = "The number of EXTERNAL mean-field averages. It is an Integer. Default is 1."
    default     = 1
    required    = false

    "--numIntAv"
    arg_type    = Int64
    help        = "The number of INTERNAL mean-field averages. It is an Integer. Default is 100."
    default     = 100
    required    = false

    "--pCCEOrder"
    arg_type    = Int64
    help        = "The order of approximation in the CCE. It is an Integer. Default is 2."
    default     = 2
    required    = false

    "--partitionSize"
    arg_type    = Int64
    help        = "The number of spins in a partition, i.e., the size of the partitions. It is an Integer. Default is 1 (this is, standard CCE)."
    default     = 1
    required    = false

    "--points"
    arg_type    = Int64
    help        = "The number of sampling points of the coherence signal. It is an Integer. Default is 25."
    default     = 25
    required    = false
    
    "--ExecuteInCluster"
    arg_type    = Bool
    help        = "True if the program will be executed in a cluster using Slurm, false otherwise (for local execution, for instance). If set to true, no need to set numWorkers, as this will be indicated in Slurm script"
    default     = false
    required    = false

    "--seed"
    arg_type    = Int64
    help        = "The seed for random number generation. It is an Integer. Set to 0 for random numbers, set otherwise for debugging. Default is 0."
    default     = 0
    required    = false

  end

  return parse_args(s; as_symbols=true)

end



function conf_driving(Species :: String, Protocol :: String)

    """
        This function configures the bath driving protocol. 

        Input parameters:
            - Species       :: The chemical species of the spin bath.
            - Protocol      :: The bath driving protocol considered.

        Output:
            - DrivenLines   :: A Dict{Int64, Float64} object where the keys are labels to the different resonance lines of the bath and the values are Bool (true if the line is driven).
            - LGDriven      :: A Dict{Int64, Float64} object where the keys are labels to the different resonance lines of the bath and the values are Bool (true if the line is LG driven).
    """

    LGDriven    = Dict(i => 0.0 for i in 1:7)
    DrivenLines = Dict(i => 0.0 for i in 1:7)

    if Protocol == "Free" return LGDriven, DrivenLines; end

    if Species == "N14"

        if Protocol == "OneToneRes"
            DrivenLines[7] = 1.0

        elseif Protocol == "TwoToneRes"
            DrivenLines[1] = DrivenLines[2] = 1.0

        elseif Protocol == "FullToneRes"
            DrivenLines[1] = DrivenLines[2] = DrivenLines[5] = DrivenLines[6] = DrivenLines[7] = 1.0

        elseif Protocol == "OneToneLG"
            LGDriven[7] = 1.0 
            DrivenLines[7] = 1.0

        elseif Protocol == "TwoToneLG"
            LGDriven[1] = LGDriven[2] = 1.0
            DrivenLines[1] = DrivenLines[2] = 1.0

        elseif Protocol == "FullToneLG"
            LGDriven[1] = LGDriven[2] = LGDriven[5] = LGDriven[6] = LGDriven[7] = 1.0
            DrivenLines[1] = DrivenLines[2] = DrivenLines[5] = DrivenLines[6] = DrivenLines[7] = 1.0

        elseif Protocol == "TwoToneHyb"
            LGDriven[2] = 1.0
            DrivenLines[1] = DrivenLines[2] = 1.0

        end

    elseif Species == "N15"

        if Protocol == "OneToneRes"
            DrivenLines[3] = 1.0

        elseif Protocol == "TwoToneRes"
            DrivenLines[3] = DrivenLines[4] = 1.0

        elseif Protocol == "FullToneRes"
            DrivenLines[1] = DrivenLines[2] = DrivenLines[3] = DrivenLines[4] = 1.0

        elseif Protocol == "OneToneLG"
            LGDriven[3] = 1.0
            DrivenLines[3] = 1.0

        elseif Protocol == "TwoToneLG"
            LGDriven[3] = LGDriven[4] = 1.0
            DrivenLines[3] = DrivenLines[4] = 1.0

        elseif Protocol == "FullToneLG"
            LGDriven[1] = LGDriven[2] = LGDriven[3] = LGDriven[4] = 1.0
            DrivenLines[1] = DrivenLines[2] = DrivenLines[3] = DrivenLines[4] = 1.0

        elseif Protocol == "TwoToneHyb"
            LGDriven[3] = 1.0
            DrivenLines[3] = DrivenLines[4] = 1.0

        end

    elseif Species == "N1415"

        if Protocol == "OneToneRes"
            DrivenLines[7] = 1.0

        elseif Protocol == "TwoToneRes"
            DrivenLines[1] = DrivenLines[2] = 1.0

        elseif Protocol == "FullToneRes"
            for i in 1:7 DrivenLines[i] = 1.0; end

        elseif Protocol == "OneToneLG"
            LGDriven[7] = 1.0
            DrivenLines[7] = 1.0

        elseif Protocol == "TwoToneLG"
            LGDriven[1] = LGDriven[2] = 1.0
            DrivenLines[1] = DrivenLines[2] = 1.0

        elseif Protocol == "TwoToneHyb"
            LGDriven[2] = 1.0
            DrivenLines[1] = DrivenLines[2] = 1.0

        elseif Protocol == "FullToneLG"
            for i in 1:7 LGDriven[i] = 1.0; DrivenLines[i] = 1.0; end
        end

    elseif Species == "e"

        if Protocol in ("OneToneRes", "TwoToneRes", "FullToneRes")
            DrivenLines[1] = 1.0
        else
            LGDriven[1] = 1.0
            DrivenLines[1] = 1.0
        end 
    end

    return LGDriven, DrivenLines
end



function summary()

    """
        This functions prints a summary of the relevant simulation parameters.
    """

    println(" ")
    println("PARAMETER SUMMARY:")
    println(" ")

    println("Number of workers:                             $(size(workers(), 1))")
    println(" ")
    
    println("Calculation under pCCE($(pCCEOrder), $(partitionSize)) approximation.")
    println(" ")

    println("Chemical species of bath spins:                $Species;")
    println("Number of bath spins in the lattice:           $NumberSpins;")
    println("Concentration of bath spin:                    [N] = $(round(concentration * 1e6, digits = 2)) ppm.")
    println(" ")

    println("Dipole radius:                                 $(round(rDipole * 1e9, digits = 2)) nm;")
    println("Mean field radius:                             $(round(rMF * 1e9, digits = 2)) nm.")
    println(" ")

    println("Number of π-pulses in the DD sequence:         $(nPulses == 0 ? "No DD sequence considered." : nPulses).")
    println(" ")

    println("Bath driving protocol:                         $(Protocol == "Free" ? "None" : Procol);")
    println("Rabi frequency of the driving:                 $(round(Int, Omega * 1e-6)) MHz;")
    println("Lines driven resonant driving:                 $(DrivenLines);")
    println("Lines driven with LG driving:                  $(LGDriven).")
    println(" ")

    println("Free evolution period duration (τ):            $(round(Int, tauMax * 1e6)) μs;")
    println("Number of sampling points:                     $points.")
    println(" ")

    println("Number of spatial averages:                    $numSpatialAv;")
    if seed != 0
      println("     Seed for random number generator is:    $seed.")
    end
    println(" ")
    
    println("MF averaging                                   $(meanField == true ? "IS" : "IS NOT") enabled.")

    if meanField
        println("Number of external averages:                   $numExtAv")
        println("Number of internal averages:                   $numIntAv")
    end
    println(" ")

    println("===================================================================================")             

end



#################################################################################################
#################################################################################################

""" Here we remove possible discrepancies between input arguments and prepare the simulation. """

if isempty(ARGS)
    include("Input.jl")
else
    input_args = parse_command_line()
    params = (; input_args...)
    (; Species, NumberSpins, concentration, nPulses, Protocol, Omega, tauMax, numSpatialAv, meanField, numExtAv, numIntAv, pCCEOrder, partitionSize, points, ExecuteInCluster, seed) = params
end

if ExecuteInCluster numWorkers = 1; end

concentration = concentration * 1e-6

rDipole                 :: Float64  = 65e-9 * (1 / (concentration * 1e6))^(1/3)
rMF                     :: Float64  = 10 * rDipole

Sparse                  :: Bool     =  partitionSize < 4 ? false : true # set true for big sparse matrices

rhoCS                   :: Matrix{Float64} = 1/2 * [1 1; 1 1] # initial state of the central spin


numSpatialAv = seed != 0 ? 1 : numSpatialAv # ensure only one spatial configuration is considered when a random seed is set.
numExtAv = meanField ? numExtAv : 1 # ensure no external averages occur when mean-field is deactivated.
numIntAv = meanField ? numIntAv : 1 # ensure no internal averages occur when mean-field is deactivated.

Omega = Protocol == "Free" ? 0.0 : Omega # ensure Rabi frequency is zero for no bath driving.

Delta                   :: Float64  = Omega / sqrt(2)
Omega_eff               :: Float64  = sqrt(Delta^2 + Omega^2)

if Omega != 0 # ensure the sampling points are not resonant with the Rabi frequency of the driving

    if Protocol in ("OneToneRes", "TwoToneRes", "FullToneRes")
        sampInteger = round(Int, Omega * tauMax / points)
        tauMax = sampInteger * points / Omega   
    else
        sampInteger = round(Int, Omega_eff * tauMax / points)
        tauMax = sampInteger * points / Omega_eff  
    end

end

LGDriven, DrivenLines = conf_driving(Species, Protocol)

Omega   *= 1e6
Delta   *= 1e6
tauMax  *= 1e-6