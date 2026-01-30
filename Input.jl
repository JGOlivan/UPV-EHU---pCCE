##############################################################################################################
######################  This file contains all the input parameters for the simulation  ######################
##############################################################################################################

Species             :: String   = "N15"     # the chemical species of the bath
NumberSpins         :: Int64    = 32       # the number of bath spins
concentration       :: Float64  = 5         # the concentration of the bath spins (in ppm)
nPulses             :: Int64    = 1         # the number of DD pulses on the central spin
Protocol            :: String   = "Free"    # the bath decoupling protocol
Omega               :: Float64  = 0         # the Rabi frequency of a single-tone bath driving (in MHz)
tauMax              :: Float64  = 30        # the free evolution period duration (in μs)
numSpatialAv        :: Int64    = 40        # the number of spatial averages (i.e. number of central spins)
meanField           :: Bool     = true      # enables mean-field averaging if true
numExtAv            :: Int64    = 1         # number of external mean-field averages 
numIntAv            :: Int64    = 40       # number of internal mean-field averages
pCCEOrder           :: Int64    = 2         # order of approximation in the CCE
partitionSize       :: Int64    = 1         # size of the partitions in pCCE (1 is equivalent to standard CCE)
points              :: Int64    = 25        # number of sampling points
parallelize         :: Bool     = true      # enable parallelization
numWorkers          :: Int64    = 4         # number of parallel workers
ExecuteInCluster    :: Bool     = false     # true if the code will execute in an HPC cluster
seed                :: Int64    = 12345         # seed for random number generation

##############################################################################################################
##############################################################################################################