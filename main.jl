include("PkgInstallation.jl")

include("Libraries.jl") # useful libraries that will only work in main.jl
include("Setup.jl")

if ExecuteInCluster 
    addprocs(SlurmManager())
elseif size(workers(), 1) == 1 && parallelize
    addprocs(numWorkers)
end

BLAS.set_num_threads(1)

@everywhere include("MyConstants.jl")
@everywhere include("MyTypes.jl")
@everywhere include("Functions.jl")

@everywhere using .MyConstants
@everywhere using .MyTypes
@everywhere using .Functions


println("----------------------------------------") 
println("----------------------------------------") 
println("TIME AT START: $(Dates.format(now(), "e dd-mm-YYYY HH:MM:SS"))")
println("----------------------------------------")             
println("----------------------------------------") 

summary() # Print a summary of the input parameters.
flush(stdout)

function main()

    # Get coherence function for given input:

    coherences = HahnSignal(pCCEOrder, partitionSize, tauMax, points, numSpatialAv, concentration, NumberSpins, seed, Species, numExtAv, numIntAv, rhoCS, rDipole, rMF, nPulses, Sparse, LGDriven, DrivenLines, meanField, Omega, Delta, ExecuteInCluster)

    coherence = mean(coherences) # average over the coherence for each spatial configuration

    # Plot coherence function:

    tauVec = LinRange(0, 2 * tauMax * 1e6, points)

    p = plot(tauVec, coherence, grid = false, framestyle = :box, legend = false, lw = 3, xlabel = L"2\tau \ (\mu s)", ylabel = L"\langle S_x \rangle", tickfont = (9, "Times New Roman"), color = :darkred, minorticks = true)

    # Save the figure to current directory

    savefig(p, "./CoherencePlot_pCCE_$(pCCEOrder)_$(partitionSize)_$(Species)_Bath_$(NumberSpins)_Spins_$(round(concentration * 1e6, digits = 2))_ppm_$(nPulses)_CS_pulses_Protocol_$(Protocol)_$(round(Int, Omega * 1e-6))_MHz_$(round(Int, tauMax * 1e6))_us_$(numSpatialAv)_sp_$(numExtAv)_ext_$(numIntAv)_int_seed_$(seed).pdf")

end

main()

println("----------------------------------------") 
println("----------------------------------------") 
println("Time at end: $(Dates.format(now(), "e dd-mm-YYYY HH:MM:SS"))")
println("----------------------------------------")             
println("----------------------------------------") 
flush(stdout)

rmprocs(workers())