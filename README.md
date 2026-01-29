# Driven Spin Bath Coherence Simulation
This is a code for simulating the coherence of a central spin surrounded by a spin bath in solid-state. The code is focused on simulating the spin coherence decay of an NV center in diamond surrounded by a bath of P1 centers. The spin bath can be driven. Some functions could be used for other purposes.  

The coherence function is calculated using the Partition Cluster Correlation Expansion method proposed in Phys. Rev. B 110, L220302 (https://doi.org/10.1103/PhysRevB.110.L220302). Their code can be found at: https://github.com/walter-hahn/pCCE-FraunhoferIAF.

This is the code used in the manuscript arXiv:2512.06948 (https://doi.org/10.48550/arXiv.2512.06948). 

## Installation
1. **Download Julia**
    You can find a guide for downloading and installing Julia at the official website: https://julialang.org
2. **Clone the repository and move to the correct directory:
    ```bash/zsh
    git clone https://github.com/JGOlivan/UPV-EHU---pCCE.git
    cd pCCE - UPV:EHU
    ```

## Usage 
The file main.jl contains a function main() where the code for performing the simulation is written. This is the part of the code you should modify. Here is a simple example of usage:
```julia
# main.jl
function main()

    # Get coherence function for given input:

    coherences = HahnSignal(pCCEOrder, partitionSize, tauMax, points, numSpatialAv, concentration, NumberSpins, seed, Species, numExtAv, numIntAv, rhoCS, rDipole, rMF, nPulses, Sparse, LGDriven, DrivenLines, meanField, Omega, Delta, ExecuteInCluster) # Returns a vector of coherence functions for different spatial configurations of the bath.

    coherence = mean(coherences) # Average over the coherence for each spatial configuration.


    # Plot coherence function:

    tauVec = LinRange(0, 2 * tauMax * 1e6, points)

    p = plot(tauVec, coherence, grid = false, framestyle = :box, legend = false, lw = 3, xlabel = L"2\tau \ (\mu s)", ylabel = L"\langle S_x \rangle", tickfont = (9, "Times New Roman"), color = :darkred, minorticks = true)


    # Save the figure to current directory

    savefig(p, "./CoherencePlot_pCCE_$(pCCEOrder)_$(partitionSize)_$(Species)_Bath_$(NumberSpins)_Spins_$(round(concentration * 1e6, digits = 2))_ppm_$(nPulses)_CS_pulses_Protocol_$(Protocol)_$(round(Int, Omega * 1e-6))_MHz_$(round(Int, tauMax * 1e6))_us_$(numSpatialAv)_sp_$(numExtAv)_ext_$(numIntAv)_int_seed_$(seed).pdf")

end
```

## Help
Every function defined in this code (in "Functions.jl") contains a small description of itself; what it is used for, what are the input arguments (and their types) and what it returns. To display this help, in the julia REPL, type '?' and then type:
```bash/zsh
include("MyConstants.jl")
include("MyTypes.jl")
include("Functions.jl")
Functions.func
```

Similarly, custom types, defined in "MyTypes.jl", can also display help by typing '?' and:
```bash/zsh
include("MyConstants.jl")
include("MyTypes.jl")
MyTypes.type
```

Note that you only need to include the files above once.

## Configuration and running the code

The code requires specific Julia packages. These will be installed automatically when running the code.

There are two ways of setting up the parameters for the simulation:
1. **From the input file:**

    ```bash/zsh
    Species             :: String   = "N15"     # the chemical species of the bath
    NumberSpins         :: Int64    = 180       # the number of bath spins
    concentration       :: Float64  = 5         # the concentration of the bath spins (in ppm)
    nPulses             :: Int64    = 1         # the number of DD pulses on the central spin
    Protocol            :: String   = "Free"    # the bath decoupling protocol
    Omega               :: Float64  = 0         # the Rabi frequency of a single-tone bath driving (in MHz)
    tauMax              :: Float64  = 30        # the free evolution period duration (in μs)
    numSpatialAv        :: Int64    = 50        # the number of spatial averages (i.e. number of central spins)
    meanField           :: Bool     = true      # enables mean-field averaging if true
    numExtAv            :: Int64    = 1         # number of external mean-field averages 
    numIntAv            :: Int64    = 100       # number of internal mean-field averages
    pCCEOrder           :: Int64    = 2         # order of approximation in the CCE
    partitionSize       :: Int64    = 1         # size of the partitions in pCCE (1 is equivalent to standard CCE)
    points              :: Int64    = 25        # number of sampling points
    parallelize         :: Bool     = true      # enable parallelization
    numWorkers          :: Int64    = 4         # number of parallel workers
    ExecuteInCluster    :: Bool     = false     # true if the code will execute in an HPC cluster
    seed                :: Int64    = 0         # seed for random number generation
    ```

    After modifying the input file, you should run:
    ```bash/zsh
    julia main.jl
    ```

2. **From the command line:**
    ```bash/zsh
    julia main.jl --arg1 val1 --arg2 val2
    ````
    The possible arguments are: --Species, --NumberSpins, --concentration, --nPulses, --Protocol, --Omega, --tauMax, --numSpatialAv, --meanField, --numExtAv, --numIntAv, --pCCEOrder, --partitionSize, --points, --ExecuteInCluster, --seed.

    For displaying help about the input arguments, type:
    ```bash/zsh
    julia main.jl --help
    ```

You shouls keep in mind that:
    1. Passing the input arguments directly from the command line takes priority over the input file. 
    2. The code can be run in parallel. To do so, start julia as:
        ```bash/zsh
        julia -p numWorkers main.jl
        ```
        Your numWorkers should not exceed the number of cores of your machine.
    3. The code can be run in an HPC (High Performance Computing) cluster that uses SLURM.