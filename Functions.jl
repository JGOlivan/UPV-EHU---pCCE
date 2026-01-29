##############################################################################################################
###############        This module contains all the functions needed for the simulation        ###############
##############################################################################################################

module Functions

using ..MyTypes, ..MyConstants
export HahnSignal

using StructArrays, SparseArrays, LinearAlgebra, Statistics, Random, Combinatorics, FastExpm, JuMP, Distributed, ProgressMeter, GLPK

"""
    createSpin(Pos :: Position{<:Real}, ID :: Int64, Species :: String)

Create a Spin object with some position in the lattice and some chemical species assigned.

    Parameters:
        - Pos         :: position vector of the spin.
        - ID          :: a unique label assigned to the spin.
        - Species     :: the chemical species of the spin.

    Returns:
        - A spin object with the specified parameters.
"""
function createSpin(Pos :: Position{<:Real}, ID :: Int64, Species :: String)


    NucSpin = Species == "N14" ? 1. : Species == "N15" ? 1/2 : 0.0
    GyroRatio = Species in ("N14", "N15", "e") ? ge : 0.0
    return Spin(Pos, ID, NucSpin, Species, 1, 0, GyroRatio, false, false) # Initialize in partition 1 and branch 0 --> this will change
end



"""
    change_basis(vec :: Position{Float64})

Change of basis from (100) to (111) diamond 

    Parameters:
        - vec        :: a position vector for a spin in the lattice
    
    Returns:
        - The vector in the rotated basis
"""
function change_basis(vec :: Position{Float64})


    theta_x = atan(sqrt(2))
    theta_z = pi / 4

    cx, sx = cos(theta_x), sin(theta_x)
    cz, sz = cos(theta_z), sin(theta_z)

    Rx = [1.0 0.0 0.0; 0.0 cx -sx; 0.0 sx cx]
    Rz = [cz -sz 0.0; sz cz 0.0; 0.0 0.0 1.0]

    return Rx * (Rz * vec)

end



""" 
    generate_position(
        iterations_x    :: Int64,
        iterations_y    :: Int64,
        iterations_z    :: Int64,
        var1            :: Vector{Int64}, 
        var2            :: Vector{Int64}, 
        var3            :: Vector{Int64}, 
        choice_list     :: Vector{Vector{Int64}},
    )

Generates a random position in the lattice 
"""
function generate_position(
    iterations_x    :: Int64,
    iterations_y    :: Int64,
    iterations_z    :: Int64,
    var1            :: Vector{Int64}, 
    var2            :: Vector{Int64}, 
    var3            :: Vector{Int64}, 
    choice_list     :: Vector{Vector{Int64}},
)

    rand1 = rand(-iterations_x:iterations_x)
    rand2 = rand(-iterations_y:iterations_y)
    rand3 = rand(-iterations_z:iterations_z)

    vec = choice_list[rand(1:length(choice_list))]

    position_vector = (vec + rand1 * var1 + rand2 * var2 + rand3 * var3) / 4

    return Position(position_vector...)
end



"""  
    lattice_generator(
        concentration       :: Float64,
        number_spins        :: Int64,
        rMF                 :: Float64,
        seed                :: Int64,
        Species             :: String
    )

Generate a list of bath spins and a list of mean-field spins

    Parameters:
        - concentration     :: the concentration of the spins in the lattice
        - number_spins      :: the number of bath spins in the lattice
        - rMF               :: the radius of mean-field lattice
        - seed              :: a seed for random number generation (0 for randomness)
        - Species           :: the chemical species of the spins in the lattice
    
    Returns:
        - A list of Spin objects containing the bath spins
        - A list of Spin objects containing the mean-field spins
"""
function lattice_generator(
    concentration       :: Float64,
    number_spins        :: Int64,
    rMF                 :: Float64,
    seed                :: Int64,
    Species             :: String
)

    # Set a random seed for reproducibility if desired

    if seed != 0
        Random.seed!(seed)
    end

    rMF /= 3.57e-10

    # Set of coordinates of Carbon atoms in the unit cell (in units of 4 times the lattice constant a):

    b1 = [0, 0, 0]
    b2 = [0, 2, 2]
    b3 = [2, 0, 2]
    b4 = [2, 2, 0]
    b5 = [3, 3, 3]
    b6 = [3, 1, 1]
    b7 = [1, 3, 1]
    b8 = [1, 1, 3]

    base_list = [b1, b2, b3, b4, b5, b6, b7, b8]

    # Set of vectors that map a position onto a contiguous unit cell (in units of 4 times the lattice constant a):
    var1 = [4, 0, 0]
    var2 = [0, 4, 0]
    var3 = [0, 0, 4]

    # We fix the number of impurities in the lattice, so its size should be bigger the less the concentration of impurities.

    size_factor = (1e-6 / concentration)^(1/3)

    # Number of unit cells generated in each direction (all three are the same as we build a cubic lattice)

    iterations_x = round(Int, 500 * size_factor)
    iterations_y = round(Int, 500 * size_factor)
    iterations_z = round(Int, 500 * size_factor)

    # Number of lattice sites generated:

    num_sites = 8 * 8 * iterations_x * iterations_y * iterations_z

    # Number of impurities generated:

    num_imp = round(Int, concentration * num_sites)

    # Preallocate the Vector of Vectors that will contain every position generated

    all_positions = Vector{Position{Float64}}(undef, num_imp)

    # Generate num_imp positions and add them to all_positions only if the position i) is not the origin (i.e. the NV) and ii) has not been previously generated:

    for n in 1:num_imp

        real_pos = generate_position(iterations_x, iterations_y, iterations_z, var1, var2, var3, base_list)

        while all(getfield(real_pos, field) == 0.0 for field in fieldnames(Position)) || any(all(real_pos .== pos) for pos in all_positions) || norm(real_pos) < 3*3.57

            real_pos = generate_position(iterations_x, iterations_y, iterations_z, var1, var2, var3, base_list)
        end

        all_positions[n] = change_basis(real_pos)
    end

    # Sort the generated positions by norm

    all_positions = sort(all_positions, by=norm)


    # Create the Spin objects for bath spins:
    bathSpins = [createSpin(pos, n, Species) for (n, pos) in enumerate(all_positions[1:number_spins])]

    max_dist = norm(all_positions[1:number_spins][end])

    # Create the Spin objects for mean-field spins:
    mfSpins = [createSpin(pos, n, Species) for (n, pos) in enumerate(all_positions) if norm(pos) <= max_dist + rMF]

    # return the bath and the mean-field

    return StructArray(bathSpins), StructArray(mfSpins)
    
end



"""
    constrainedKMeans(
        X            :: Vector{Position{Float64}}, 
        K            :: Int64,
        tau_max      :: Vector{Int64},
        tau_min      :: Vector{Int64},
        maxiter      :: Int64
    )

This function partitions a dataset applying a constrained KMeans algorithm.

        Parameters:
            - X         :: the dataset (a list of positions)
            - K         :: the number of partitions desired
            - tau_max   :: the maximum size for each partition
            - tau_min   :: the minimum size for each partition
            - maxiter   :: the maximum number of iterations of the optimization algorithm
        
        Returns:
            - A list of labels for each data point corresponding to the assigned partition
"""
function constrainedKMeans(
    X            :: Vector{Position{Float64}}, 
    K            :: Int64,
    tau_max      :: Vector{Int64},
    tau_min      :: Vector{Int64},
    maxiter      :: Int64
    )

    # Get dataset dimensions.
    n = length(X) # number of spins in the bath 
    D = length(Position) # dimension of the space
    #n, D = size(X)

    # Randomly assign each data point to a cluster.
    labels = rand(1:K, n)

    centroids = zeros(K, D)

    # Initialize centroids.
    for k in 1:K

        #println("Centroid number: ", k)

        pts = X[labels .== k]

        if size(pts, 1) > 0
            centroids[k, :] = vec(mean(pts))
        else # if there is no data point in the cluster, take a random point to be the centroid
            centroids[k, :] = X[rand(1:n)]
        end

        # If the computed centroid is invalid, choose a random point.
        if any(isnan.(centroids[k, :])) || all(centroids[k, :] .== 0)
            centroids[k, :] = X[rand(1:n)]
        end

    end

    iter = 1

    while iter < maxiter

        if iter % 10 == 0
            println("KMeans iteration number: ", iter)
        end

        # Compute squared Euclidean distances
        # We obtain an n×K matrix.

        distmat = [sum((X[i] .- centroids[k, :]).^2) for i in 1:n, k in 1:K]

        # Flatten to a vector of length n*K.
        objcoef = reshape(distmat, n*K)

        # Set up the MILP to assign each point to one cluster.
        model = Model(GLPK.Optimizer)
        set_silent(model)  # Suppress solver output

        # Create binary variables t[i,k] for i=1:n and k=1:K.
        @variable(model, t[1:n, 1:K], Bin)

        # Each data point is assigned to exactly one cluster.
        @constraint(model, [i=1:n], sum(t[i, k] for k in 1:K) == 1)

        # Enforce maximum cluster size: for each cluster k, sum_{i} t[i,k] ≤ tau_max[k]
        @constraint(model, [k=1:K], sum(t[i, k] for i in 1:n) ≤ tau_max[k])

        # Enforce minimum cluster size: for each cluster k, sum_{i} t[i,k] ≥ tau_min[k]
        @constraint(model, [k=1:K], sum(t[i, k] for i in 1:n) ≥ tau_min[k])

        # Set objective: minimize the total squared distance.
        @objective(model, Min, sum(objcoef[(k-1)*n + i] * t[i, k] for i in 1:n, k in 1:K))

        optimize!(model)

        # Retrieve the solution. 
        T_sol = value.(t)

        # For each data point, assign the cluster corresponding to the variable set to 1.
        labelsNew = [argmax(T_sol[i, :]) for i in 1:n]

        # Compute new centroids.
        centroidsNew = zeros(K, D)

        for k in 1:K
            pts = X[labelsNew .== k]

            if size(pts, 1) > 0
                centroidsNew[k, :] = vec(mean(pts))
            else
                centroidsNew[k, :] = X[rand(1:n)]
            end
        end

        # Check for convergence

        if norm(centroids - centroidsNew) < 1e-6
            labels = labelsNew
            centroids = centroidsNew
            break
        else
            centroids = centroidsNew
            labels = labelsNew
            iter += 1
        end

    end

    return labels

end



"""
    get_partitions!(
        partitionSize       :: Int64,
        bathSpins           :: StructArray
    )

Partitions a bath of spins into groups of strongly interacting spins using a constrained KMeans algorithm.

    Parameters:
        - partition Size    :: the size of the partitions
        - bathSpins         :: the list of bath spins

    Returns:
        - the list of bath spins with the assigned partition
"""
function get_partitions!(
    partitionSize       :: Int64,
    bathSpins           :: StructArray
)

    numP1 = bathSpins.ID[end]

    numPartitions = partitionSize == 1 ? numP1 : mod(numP1, partitionSize) == 0 ? div(numP1, partitionSize) : div(numP1, partitionSize) + 1

    tauMax = fill(partitionSize, numPartitions)
    tauMin = ones(Int64, numPartitions)
    maxiter = 300

    labels = constrainedKMeans(bathSpins.Pos, numPartitions, tauMax, tauMin, maxiter)

    for i in 1:numPartitions bathSpins.Partition[labels .== i] .= i; end

    return bathSpins
end



"""
    hyperfineNV(r :: Position{<:Real}, GyroRatio :: Float64)

Computes the interaction vector betweem an NV center (or a simple e spin) and another spin (nuclear or electron).

    Parameters:
        - r             :: the position of the spin with respect to the NV center
        - GyroRatio     :: the gyromagnetic ratio of the spin

    Returns:
        - the components of the interaction vector Ax, Ay, Az.
"""
function hyperfineNV(r :: Position{<:Real}, GyroRatio :: Float64)

    G = ge * GyroRatio * mu0 * hbar / 2.0

    Ax = G / (norm(r)^3) * (-3 * (r.x * r.z) / norm(r)^2)
    Ay = G / (norm(r)^3) * (-3 * (r.y * r.z) / norm(r)^2)
    Az = G / (norm(r)^3) * (1 - 3 * (r.z * r.z) / norm(r)^2)
    
    return Ax, Ay, Az

end



"""
    DipoleDipole(r :: Position{Float64}, GyroRatio1 :: Float64, GyroRatio2 :: Float64)

Computes the interaction strength between two bath spins with similar Larmor frequency.

    Parameters:
        - r             :: the relative vector connecting the two spins
        - GyroRatio1    :: the gyromagnetic ratio of the first spin.
        - GyroRatio2    :: the gyromagnetic ratio of the second spin.

    Returs:
        - the interaction coupling C.
"""
function DipoleDipole(r :: Position{Float64}, GyroRatio1 :: Float64, GyroRatio2 :: Float64)


    G = GyroRatio1 * GyroRatio2 * mu0 * hbar / 2.0

    return G / norm(r)^3 * (1 - 3 * (r.z * r.z) / norm(r)^2)
end



""" 
    getSpinOps(N :: Int64, Sparse :: Bool)
Obtain NV and P1 spin operators for a (1 + N) spin system
"""
function getSpinOps(N :: Int64, Sparse :: Bool)

    sigmax :: Matrix{ComplexF64} = [0 1; 1 0] # these are now defined in constants
    sigmay :: Matrix{ComplexF64} = [0 -im; im 0]
    sigmaz :: Matrix{ComplexF64} = [1 0; 0 -1]

    if Sparse
        sigmax, sigmay, sigmaz = sparse(sigmax), sparse(sigmay), sparse(sigmaz)
    end

    identityMat = Sparse ? sparse(1.0I(2^N)) : 1.0I(2^N)

    # NV Electron Spin operators

    Sx = kron(sigmax, identityMat)
    Sy = kron(sigmay, identityMat)
    Sz = kron(sigmaz, identityMat)

    opType = Sparse ? SparseMatrixCSC{ComplexF64, Int64} : Matrix{ComplexF64}

    # P1 Electron Spin operators

    Jx = Vector{opType}(undef, N)
    Jy = similar(Jx)
    Jz = similar(Jx)

    NV = Sparse ? sparse(1.0I(2)) : 1.0I(2)

    for i in 1:N
        sL = 2^(i - 1)
        sR = 2^(N - i)

        L = Sparse ? sparse(1.0I(sL)) : 1.0I(sL)
        R = Sparse ? sparse(1.0I(sR)) : 1.0I(sR)

        Jx[i] = kron(kron(kron(NV, L), 0.5 * sigmax), R)
        Jy[i] = kron(kron(kron(NV, L), 0.5 * sigmay), R)
        Jz[i] = kron(kron(kron(NV, L), 0.5 * sigmaz), R)
    end

    return spinOperators(Sx, Sy, Sz, Jx, Jy, Jz)
end



""" 
    GetBranches!(bathSpins :: StructArray)

Assigns a resonance line to the bath spins.

    Parameters:
        - bathSpins     :: the list of bath spins

    Returns:
        - the list of bath spins with assigned branches

    The following is a schematic depiction of the P1 bath resonance structure:

    14N P1 bath:        6 2 7 1 5
                        | | | | |
                          | | |
                          | | |
                            |

    15N P1 bath:        2 4 3 1
                        | | | |
                          | |
                          | |

    14N + 15N P1 bath:  6 2 4 7 3 1 5
                        | | | | | | |
                          | | | | |
                          | | | | |
                          |   |   |
"""
function GetBranches!(bathSpins :: StructArray)

    if bathSpins.Species[1] == "e" # Electrons: 1 branch

        for i in eachindex(bathSpins)
            bathSpins.Branch[i] = 1
        end

    elseif bathSpins.Species[1] == "N14" # Nitrogen 14: 5 branches

        for i in eachindex(bathSpins)
            r = rand()
            bathSpins.Branch[i] = r <= 1/12 ? 5 : r <= 2/12 ? 6 : r<= 5/12 ? 1 : r<= 8/12 ? 2 : 7
        end

    elseif bathSpins.Species[1] == "N15" # Nitrogen 15: 4 branches

        for i in eachindex(bathSpins)
            r = rand()
            bathSpins.Branch[i] = r <= 1/8 ? 1 : r <= 2/8 ? 2 : r <= 5/8 ? 3 : 4
        end

    elseif bathSpins.Species[1] == "N1415" # Nitrogen mixture: 7 branches

        if rand() <= 0.996

            for i in eachindex(bathSpins)
                r = rand()
                bathSpins.Branch[i] = r <= 1/12 ? 5 : r <= 2/12 ? 6 : r <= 5/12 ? 1 : r <= 8/12 ? 2 : 7
            end

        else

            for i in eachindex(bathSpins)
                r = rand()
                bathSpins.Branch[i] = r <= 1/8 ? 1 : r <= 2/8 ? 2 : r <= 5/8 ? 3 : 4
            end
        end

    end

    return bathSpins
end



"""
    DrivenSpins!(
        bathSpins       :: StructArray,
        DrivenLines     :: Dict{Int64, Float64},
        LGDriven        :: Dict{Int64, Float64}
    )

Assigns a true or false label to driven/not driven spins based on the asssigned branch.

    Parameters:
        - bathSpins     :: the list of bath spins.
        - DrivenLines   :: a dictionary that indicates which lines are driven.
        - LDGriven      :: a dictionary that indicates which lines are LG driven.

    Returns:
        - the list of bath spins with the corresponding labels
"""
function DrivenSpins!(
    bathSpins       :: StructArray,
    DrivenLines     :: Dict{Int64, Float64},
    LGDriven        :: Dict{Int64, Float64}
)

    for i in eachindex(bathSpins)
        bathSpins.Driven[i] = DrivenLines[bathSpins.Branch[i]] == 1.0 ? true : false
        bathSpins.LGDriven[i] = LGDriven[bathSpins.Branch[i]] == 1.0 ? true : false
    end

    return bathSpins
end



"""
    getHamiltonian(
        indices             :: Vector{Int64},
        D                   :: Driving{Float64, String},
        rDipole             :: Float64,
        Sparse              :: Bool,
        spinOps             :: spinOperators,
        bathSpins           :: StructArray
    )

This function obtains the Hamiltonian operator for some partition(s) of bath spins.

    Parameters:
        - indices       :: A list of indices corresponding to some partition(s).
        - D             :: A Driving type describing a single bath driving tone.
        - rDipole       :: A cutoff distance for including interactions between bath spins.
        - Sparse        :: True for working with sparse matrices.
        - spinOps       :: The spin operators for every spin.
        - bathSpins     :: The list of all bath spins.

    Returns:
        - The P1 Hamiltonian.
        - The P1-NV interaction Hamiltonian.
        - The P1-P1 interaction Hamiltonian within a partition.
        - The driving Hamiltonian along A axis (Resonant and LG driving). 
        - The driving Hamiltonian along B axis (FSLG and LG4 driving).
        - The P1-P1 interaction between partitions (only for more than one partition).
"""
function getHamiltonian(
    indices             :: Vector{Int64},
    D                   :: Driving{Float64, String},
    rDipole             :: Float64,
    Sparse              :: Bool,
    spinOps             :: spinOperators,
    bathSpins           :: StructArray
)

    alpha = D.alpha
    Omega = D.Omega
    Delta = D.Delta
    LG = D.LG

    Sz, Jx, Jy, Jz = spinOps.Sz, spinOps.Jx, spinOps.Jy, spinOps.Jz

    one = Sparse ? sparse(I, size(Sz)...) : Matrix{ComplexF64}(I, size(Sz)...)

    HP1   = spzeros(ComplexF64, size(Sz)...)
    HNV   = spzeros(ComplexF64, size(Sz)...)
    HP1P1 = spzeros(ComplexF64, size(Sz)...)
    H12   = spzeros(ComplexF64, size(Sz)...)
    HA    = spzeros(ComplexF64, size(Sz)...)
    HB    = spzeros(ComplexF64, size(Sz)...)

    startIdx = 0

    for idx in sort(indices)

        mask = bathSpins.Partition .== idx
        partitionIDs = bathSpins.ID[mask]

        partitionSize = size(partitionIDs, 1)

        HP1Partition   = spzeros(ComplexF64, size(Sz)...)
        HNVPartition   = spzeros(ComplexF64, size(Sz)...)
        HP1P1Partition = spzeros(ComplexF64, size(Sz)...)
        H12Partition   = spzeros(ComplexF64, size(Sz)...)
        HAPartition    = spzeros(ComplexF64, size(Sz)...)
        HBPartition    = spzeros(ComplexF64, size(Sz)...)

        for (s, id) in enumerate(partitionIDs)

            JXs = Jx[startIdx + s]
            JYs = Jy[startIdx + s]
            JZs = Jz[startIdx + s]

            HP1Partition    += bathSpins.LGDriven[id] * Delta * JZs
            HAPartition     += bathSpins.Driven[id] * Omega * (sin(alpha) * JXs + cos(alpha) * JYs)

            if LG == "LG4" HBPartition     += bathSpins.Driven[id] * Omega * (sin(-alpha) * JXs + cos(-alpha) * JYs); end

            _, _, Az = hyperfineNV(bathSpins.Pos[id] * a, bathSpins.GyroRatio[id])

            HNVPartition += Az * (Sz + one) / 2 * JZs

            for r in 1:(s - 1) # Intra partition interactions

                id2 = partitionIDs[r]

                JXr = Jx[startIdx + r]
                JYr = Jy[startIdx + r]
                JZr = Jz[startIdx + r]

                r12 = bathSpins.Pos[id] - bathSpins.Pos[id2]
                C12 = DipoleDipole(r12 * a, bathSpins.GyroRatio[id], bathSpins.GyroRatio[id2]) # Positions are given in lattice constant units. This could be different for lattices other than diamond

                if norm(r12 * a) < (1.5 * 1.54e-10) 
                    continue # Two adjacent P1 centres (dist < bond length) will create a different kind of defect
                end

                HP1P1Partition += C12 * (JZs * JZr - (bathSpins.Branch[id] == bathSpins.Branch[id2] ? 0.5 : 0.0) * (JXs * JXr + JYs * JYr))

            end

            # Include inter-partition interactions

            startIdx2 = partitionSize

            for idx2 in filter(x -> x != idx && x > idx, indices)

                mask2 = bathSpins.Partition .== idx2

                partitionIDs2 = bathSpins.ID[mask2]
                partitionSize2 = size(partitionIDs2, 1)

                for (u, id2) in enumerate(partitionIDs2)

                    JXu = Jx[startIdx2 + u]
                    JYu = Jy[startIdx2 + u]
                    JZu = Jz[startIdx2 + u]

                    r12 = bathSpins.Pos[id] - bathSpins.Pos[id2]

                    if norm(r12 * a) < (1.5 * 1.54e-10) || norm(r12 * a) > rDipole continue; end # i) Two adjacent P1 centres (dist < bond length) will create a different kind of defect. ii) Interaction between two P1 centres separated further than the dipole radius are considered to be suppressed.

                    C12 = DipoleDipole(r12 * a, bathSpins.GyroRatio[id], bathSpins.GyroRatio[id2])

                    H12Partition += C12 * (JZs * JZu - (bathSpins.Branch[id] == bathSpins.Branch[id2] ? 0.5 : 0.0) * (JXs * JXu + JYs * JYu))

                end
                startIdx2 += partitionSize2
            end
        end

        HP1     += HP1Partition
        HNV     += HNVPartition
        HP1P1   += HP1P1Partition
        HA      += HAPartition
        HB      += HBPartition
        H12     += H12Partition

        startIdx += partitionSize
    end

    return HP1, HNV, HP1P1, HA, HB, H12

end



"""
    getMFHamiltonian(
        bathSpins               :: StructArray,
        mfSpins                 :: StructArray,
        indices                 :: Vector{Int64},
        mfStates                :: Vector{Float64},
        totalPositions          :: Vector{Position{Float64}},
        spinOps                 :: spinOperators,
        D                       :: Driving{Float64, String}  
    )

This function obtains the mean-field part of the Hamiltonian for some partition(s) of bath spins.

    Parameters:
        - bathSpins         :: The list of all bath spins.
        - mfSpins           :: The list of all mean-field spins.
        - indices           :: A list of indices corresponding to some partition(s).
        - mfStates          :: A list of random states for mean-field spins.
        - totalPositions    :: A list of positions of the spin in the partition(s).
        - spinOps           :: The spin operators for all the spins.
        - D                 :: A Driving object describing a single tone of bath driving.
        
    Returns:
        - The mean-field part of the Hamiltonian.
"""
function getMFHamiltonian(
    bathSpins               :: StructArray,
    mfSpins                 :: StructArray,
    indices                 :: Vector{Int64},
    mfStates                :: Vector{Float64},
    totalPositions          :: Vector{Position{Float64}},
    spinOps                 :: spinOperators,
    D                       :: Driving{Float64, String}  
)

    alpha = D.alpha

    Sz, Jx, Jy, Jz = spinOps.Sz, spinOps.Jx, spinOps.Jy, spinOps.Jz

    HMF = spzeros(ComplexF64, size(Sz)...)

    startIdx = 0

    for idx in sort(indices)
        
        mask = bathSpins.Partition .== idx
        
        partitionIDs = bathSpins.ID[mask]
        partitionSize = size(partitionIDs, 1)
        
        HMFPartition = spzeros(ComplexF64, size(Sz)...)

        for (s, id) in enumerate(partitionIDs)

            JXs = Jx[startIdx + s]
            JYs = Jy[startIdx + s]
            JZs = Jz[startIdx + s]

            meanFieldZ = 0
            meanFieldAlpha = 0
            meanFieldP = 0

            spinDriven          = bathSpins.Driven[id]
            spinLGDriven        = bathSpins.LGDriven[id]
            spinBranch          = bathSpins.Branch[id]
            spinPosition        = bathSpins.Pos[id]
            
            for (p, mfspin) in enumerate(mfSpins)

                if mfspin.Pos in totalPositions || (spinBranch == mfspin.Branch && spinLGDriven * mfspin.LGDriven) # Explicitly included spins and double magic pairs do not have a mean-field
                    continue
                elseif (spinLGDriven == mfspin.LGDriven && spinDriven != mfspin.Driven)
                    r = spinPosition - mfspin.Pos
                    if norm(r * a) < (1.5 * 1.54e-10) continue; end
                    C = DipoleDipole(r * a, bathSpins.GyroRatio[s], mfspin.GyroRatio)

                    if spinDriven meanFieldZ += C * mfStates[p]; end

                    continue
                end

                r = spinPosition - mfspin.Pos

                if norm(r * a) < (1.5 * 1.54e-10) continue; end # P1s that are closer than the diamond bond length are assumed to form a different defect and skipped

                C = DipoleDipole(r * a, bathSpins.GyroRatio[s], mfspin.GyroRatio) # we dont care about which ratio to take because both are the same in this case

                if mfspin.Driven * spinDriven && spinLGDriven == mfspin.LGDriven # both resonant, same or different branch, or both LG different branch

                    if spinLGDriven == true # both LG (only in different branchs!)
                        meanFieldP += C * 1/3 * mfStates[p] # it wont enter here because we do not consider this case for now
                    elseif spinBranch == mfspin.Branch # both resonant, same branch
                        meanFieldAlpha += -C * mfStates[p] * 1/2
                    else # both resonant, different branch
                        continue
                    end

                elseif spinLGDriven && !(mfspin.Driven) # bath spin magic, the other none 
                    meanFieldP += C * mfStates[p] * 1/sqrt(3) # cos(theta) = 1 / √3

                elseif mfspin.LGDriven && !(spinDriven) # bath spin none, the other magic
                    meanFieldZ += C * mfStates[p] * 1/sqrt(3)

                elseif spinDriven * mfspin.Driven && (spinLGDriven || mfspin.LGDriven) # one magic, the other resonant
                    continue

                else # other cases (none driven) normal mean field
                    meanFieldZ += C * mfStates[p] # mfstates is +- 1/2

                end
            end

            HMFPartition += meanFieldAlpha * (sin(alpha) * JXs + cos(alpha) * JYs)
            HMFPartition += meanFieldZ * JZs
            HMFPartition += meanFieldP * (1/sqrt(3) * JZs + sqrt(2/3) * (sin(alpha) * JXs + cos(alpha) * JYs))

        end

        HMF += HMFPartition
        startIdx += partitionSize
    end

    return HMF
end



"""
    GetIndices(order :: Int64, N :: Int64)

Gets the indices corresponding to each term in the CCE at some order. 

    Example for a 4 spin system (N=4):
        - order = 1: returns: [[1], [2], [3], [4]]
        - order = 2: returns: [[2, 1], [3, 1], [4, 1], [3, 2], [4, 2], [4, 3]]
        - order = 3: returns: [[3, 2, 1], [4, 2, 1], [4, 3, 1], [4, 3, 2]]
        - order = 4: returns: [[4, 3, 2, 1]]
"""
function GetIndices(order :: Int64, N :: Int64)


    ranges = [range(order - k + 1, N + 1 - k) for k in 1:order]
    allIndices = collect(Iterators.product(ranges...))
    indices = [collect(tuple) for tuple in allIndices if all(tuple[i] < tuple[i-1] for i in 2:order)]
    return indices
end


"""
    IntAvSignal(
        mfStates            :: Vector{Float64},
        bathSpins           :: StructArray,      
        mfSpins             :: StructArray,
        index               :: Vector{Int64},
        totalPositions      :: Vector{Position{Float64}},
        spinOps             :: spinOperators,
        D                   :: Driving{Float64, String},
        H0                  :: Matrix{ComplexF64},
        timeVec             :: AbstractVector,
        ux                  :: Matrix{ComplexF64},
        nPulse              :: Int64,
        rho0                :: Matrix{Float64},
        meanField           :: Bool
    )

This function computes a coherence term in the CCE for one sample of the mean-field.

    Parameters:
        - mfStates          :: A vector of random states for the mean-field spins.
        - bathSpins         :: The list of bath spins.
        - mfSpins           :: The list of mean-field spins.
        - index             :: The indices corresponding to the considered partition(s).
        - totalPositions    :: The positions of all the spins in the considered partition(s).
        - spinOps           :: The spin operators of all the spins.
        - D                 :: A Driving type describing a single tone bath driving.
        - H0                :: The Hamiltonian operator for the considered partition(s).
        - timeVec           :: The vector containing the sampled of tau vectors.
        - ux                :: The pi-pulse operator for the central spin.
        - nPulse            :: The number of pulses applied to the central spin.
        - rho0              :: The initial state of the central spin.
        - meanField         :: True if a mean field is to be considered.

    Returns:
        - A coherence term in the CCE for a single sample of the mean-field.
"""
function IntAvSignal(
    mfStates            :: Vector{Float64},
    bathSpins           :: StructArray,      
    mfSpins             :: StructArray,
    index               :: Vector{Int64},
    totalPositions      :: Vector{Position{Float64}},
    spinOps             :: spinOperators,
    D                   :: Driving{Float64, String},
    H0                  :: Matrix{ComplexF64},
    timeVec             :: AbstractVector,
    ux                  :: Matrix{ComplexF64},
    nPulse              :: Int64,
    rho0                :: Matrix{Float64},
    meanField           :: Bool
    )

    for p in eachindex(mfStates) mfStates[p] = rand() > 1/2 ? 1/2 : -1/2; end

    meanField ? HMF = getMFHamiltonian(bathSpins, mfSpins, index, mfStates, totalPositions, spinOps, D) : HMF = zeros(size(H)...)

    H = H0 + HMF

    U0 = sparse(fastExpm(-im * 2*pi * H * timeVec[2]; threshold = 1e-6)) 

    U = 1.0I(size(U0, 1))
    sig = zeros(length(timeVec))

    for i in eachindex(timeVec)
        nPulse == 0 ? U1 = U * U : U1 = (U * ux * U)^nPulse
        sig[i] += real(tr(U1 * rho0 * U1' * spinOps.Sx))
        U = U0 * U
    end

    return sig

end


"""
    pCCECalculation(
        pCCEOrder           :: Int64,
        rDipole             :: Float64,
        numExtAv            :: Int64,
        numIntAv            :: Int64,
        bathSpins           :: StructArray,
        mfSpins             :: StructArray,
        tauMax              :: Float64,
        points              :: Int64,
        D                   :: Driving{Float64, String},
        nPulse              :: Int64,
        rhoNV               :: Matrix{Float64},
        Sparse              :: Bool,
        meanField           :: Bool
    )

This function computes the coherence function for a single bath configuration.

    Parameters:
        - pCCEOrder         :: The order of approximation in the CCE.
        - rDipole           :: The cutoff distance for considering interactions between bath spins.
        - numExtAv          :: The number of external mean-field averages.        
        - numIntAv          :: The number of internal mean-field averages.  
        - bathSpins         :: The list of bath spins.
        - mfSpins           :: The list of mean-field spins.
        - tauMax            :: The free evolution period duration.
        - points            :: The number of sampling points for the coherence function.
        - D                 :: A Driving type describing a single bath driving tone.
        - nPulse            :: The number of pi-pulses applied to the central spin.
        - rhoNV             :: The initial state of the central spin.
        - Sparse            :: True for using sparse matrices.
        - meanField         :: True for including a mean-field.  
        
    Returns:
        - The coherence function for a single bath configuration.
"""
function pCCECalculation(
    pCCEOrder           :: Int64,
    rDipole             :: Float64,
    numExtAv            :: Int64,
    numIntAv            :: Int64,
    bathSpins           :: StructArray,
    mfSpins             :: StructArray,
    tauMax              :: Float64,
    points              :: Int64,
    D                   :: Driving{Float64, String},
    nPulse              :: Int64,
    rhoNV               :: Matrix{Float64},
    Sparse              :: Bool,
    meanField           :: Bool
)

    averagedSignal = zeros(points)
    timeVec = LinRange(0, tauMax, points)

    numPartitions = maximum(bathSpins.Partition)

    identityFunc(n) = Sparse ? sparse(I(2*2^n)) : 1.0I(2*2^n)
    rho0P1(n) = Sparse ? sparse(1/2^n * I(2^n)) : (1/2^n) * 1.0I(2^n)

    mfStates = zeros(mfSpins.ID[end])

    signal = ones(Float64, points)

    for _ in 1:numExtAv

        fill!(signal, 1.0)

        prevSignal = ones(points, numPartitions)
        prevLabels = Vector{Vector{Int64}}()

        for order in 1:pCCEOrder

            indices = GetIndices(order, numPartitions)
            numContributions = length(indices)

            Labels = [zeros(Int64, order) for _ in 1:numContributions]
            Signal = [zeros(Float64, points) for _ in 1:numContributions]

            for (counter, index) in enumerate(indices)

                mask = in(index).(bathSpins.Partition)
                totalSize = length(bathSpins.Partition[mask])
                totalPositions = bathSpins.Pos[mask]

                spinOps = getSpinOps(totalSize, Sparse)
                ux = -im * spinOps.Sx
                rho0 = kron(rhoNV, rho0P1(totalSize))

                HP1, HNV, HP1P1, HA, HB, H12 = getHamiltonian(index, D, rDipole, Sparse, spinOps, bathSpins)

                H0 = HP1 + HNV + HP1P1 + H12 + HA + HB

                sig = zeros(points)

                for _ in 1:numIntAv
                    sig .+= IntAvSignal(mfStates, bathSpins, mfSpins, index, totalPositions, spinOps, D, H0, timeVec, ux, nPulse, rho0, meanField)
                end

                sig /= numIntAv

                if order > 1

                    prevSignalIndices = findall(lab -> lab in combinations(index, order-1), prevLabels)

                    denominator = prod(prevSignal[:, prevSignalIndices], dims = 2)

                    Signal[counter] .= sig ./ vec(denominator)

                else
                    Signal[counter] .= sig
                end

                Labels[counter] = index
            end

            if !isempty(Signal)
                newSignals = reduce(hcat, Signal)
                signal .*= prod(newSignals, dims = 2)
            end

            prevSignal = newSignals
            prevLabels = Labels
        end

        averagedSignal .+= signal
    end

    return averagedSignal ./ numExtAv
end


"""
    Compute(
        concentration       :: Float64,
        numP1               :: Int64,
        rMF                 :: Float64,
        seed                :: Int64,
        partitionSize       :: Int64,
        Species             :: String,
        pCCEOrder           :: Int64,
        rDipole             :: Float64,
        numExtAv            :: Int64,
        numIntAv            :: Int64,
        tauMax              :: Float64,
        points              :: Int64,
        Omega               :: Float64,
        Delta               :: Float64,
        nPulse              :: Int64,
        LGDriven            :: Dict{Int64, Float64},
        DrivenLines         :: Dict{Int64, Float64},
        rhoNV               :: Matrix{Float64},
        Sparse              :: Bool,
        meanField           :: Bool
    )

This function generates a bath configuration and computes the coherence function.

    Parameters:
        - concentration     :: The concentration of the bath spins.
        - numP1             :: The number of bath spins.
        - rMF               :: The radius of the mean-field sphere.
        - seed              :: A seed for random number generation (set different to 0 for debugging).
        - partitionSize     :: The size of the partitions in pCCE.
        - Species           :: The chemical species of the bath spins.
        - pCCEOrder         :: The order of approximation in the CCE.
        - rDipole           :: The cutoff distance for considering interactions between bath spins.
        - numExtAv          :: The number of external mean-field averages.
        - numIntAv          :: The number of internal mean-field averages.
        - tauMax            :: The free evolution period duration.
        - points            :: The number of sampling points for the coherence function.
        - Omega             :: The Rabi frequency of the bath driving.
        - Delta             :: The detuning of the bath driving.
        - nPulse            :: The number of pi-pulses applied to the central spin.
        - LGDriven          :: A dictionary indicating the resonance lines LG driven.
        - DrivenLines       :: A dictionary indicating the resonance lines driven.
        - rhoNV             :: The initial state of the central spin.
        - Sparse            :: True for sparse matrices.
        - meanField         :: True for including a mean-field.

    Returns:
        - The coherence function for a single spatial configuration of the bath.
 """
function Compute(
    concentration       :: Float64,
    numP1               :: Int64,
    rMF                 :: Float64,
    seed                :: Int64,
    partitionSize       :: Int64,
    Species             :: String,
    pCCEOrder           :: Int64,
    rDipole             :: Float64,
    numExtAv            :: Int64,
    numIntAv            :: Int64,
    tauMax              :: Float64,
    points              :: Int64,
    Omega               :: Float64,
    Delta               :: Float64,
    nPulse              :: Int64,
    LGDriven            :: Dict{Int64, Float64},
    DrivenLines         :: Dict{Int64, Float64},
    rhoNV               :: Matrix{Float64},
    Sparse              :: Bool,
    meanField           :: Bool
)

    bathSpins, mfSpins  = lattice_generator(concentration, numP1, rMF, seed, Species)
    bathSpins           = get_partitions!(partitionSize, bathSpins)

    bathSpins           = GetBranches!(bathSpins)
    bathSpins           = DrivenSpins!(bathSpins, DrivenLines, LGDriven)

    mfSpins             = GetBranches!(mfSpins)

    for id in bathSpins.ID mfSpins.Branch[id] = bathSpins.Branch[id]; end

    mfSpins             = DrivenSpins!(mfSpins, DrivenLines, LGDriven)

    D                   = Driving(Omega, Delta, 0.0, "LG")

    signal              = pCCECalculation(pCCEOrder, rDipole, numExtAv, numIntAv, bathSpins, mfSpins, tauMax, points, D, nPulse, rhoNV, Sparse, meanField)

    return signal
end

"""
    HahnSignal(
        pCCEOrder           :: Int64,
        partitionSize       :: Int64,
        tauMax              :: Float64,
        points              :: Int64,
        K                   :: Int64,
        concentration       :: Float64,
        numP1               :: Int64,
        seed                :: Int64,
        Species             :: String,
        numExtAv            :: Int64,
        numIntAv            :: Int64,
        rhoNV               :: Matrix{Float64},
        rDipole             :: Float64,
        rMF                 :: Float64,
        nPulse              :: Int64,
        Sparse              :: Bool,
        LGDriven            :: Dict{Int64, Float64},
        DrivenLines         :: Dict{Int64, Float64},
        meanField           :: Bool,
        Omega               :: Float64,
        Delta               :: Float64,
        ExecuteInCluster    :: Bool
    )

This function computes the coherence function for a number of different spatial configurations of the spin bath. This is the main function of the program.

    Parameters:
        - pCCEOrder         :: The order of approximation in the CCE.
        - partitionSize     :: The size of the partitions in pCCE.
        - tauMax            :: The free evolution period duration.
        - points            :: The number of sampling points for the coherence function.
        - K                 :: The number of spatial configurations considered.
        - concentration     :: The concentration of the bath spins.
        - numP1             :: The number of bath spins.
        - seed              :: A seed for random number generation (set different to 0 for debugging).
        - Species           :: The chemical species of the bath spins.
        - numExtAv          :: The number of external mean-field averages.
        - numIntAv          :: The number of internal mean-field averages.
        - rhoNV             :: The initial state of the central spin.
        - rDipole           :: The cutoff distance for considering interactions between bath spins.
        - rMF               :: The radius of the mean-field sphere.
        - nPulse            :: The number of pi-pulses applied to the central spin.
        - Sparse            :: True for sparse matrices.
        - LGDriven          :: A dictionary indicating the resonance lines LG driven.
        - DrivenLines       :: A dictionary indicating the resonance lines driven.
        - meanField         :: True for including a mean-field.
        - Omega             :: The Rabi frequency of the bath driving.
        - Delta             :: The detuning of the bath driving.
        - ExecuteInCluster  :: True if the code is to be executed in an HPC cluster.

    Returns:
        - The coherence function for a number of spatial configurations of the bath.
 """
function HahnSignal(rhoNV :: Matrix{Float64}, rDipole :: Float64, rMF :: Float64;
    pCCEOrder           :: Int64                    = 2,
    partitionSize       :: Int64                    = 1,
    tauMax              :: Float64                  = 30e-6,
    points              :: Int64                    = 25,
    K                   :: Int64                    = 100,
    concentration       :: Float64                  = 5e-6,
    numP1               :: Int64                    = 180,
    seed                :: Int64                    = 0,
    Species             :: String                   = "N15",
    numExtAv            :: Int64                    = 1,
    numIntAv            :: Int64                    = 100,
    nPulse              :: Int64                    = 1,
    Sparse              :: Bool                     = false,
    LGDriven            :: Dict{Int64, Float64}     = Dict(i => 0.0 for i in 1:7),
    DrivenLines         :: Dict{Int64, Float64}     = Dict(i => 0.0 for i in 1:7),
    meanField           :: Bool                     = true,
    Omega               :: Float64                  = 0.0,
    Delta               :: Float64                  = 0.0,
    ExecuteInCluster    :: Bool                     = false
)

    if size(workers(), 1) > 1
        if ExecuteInCluster
            signal = pmap(k -> Compute(concentration, numP1, rMF, seed, partitionSize, Species, pCCEOrder, rDipole, numExtAv, numIntAv, tauMax, points, Omega, Delta, nPulse, LGDriven, DrivenLines, rhoNV, Sparse, meanField), 1:K)
        else
            signal = @showprogress pmap(k -> Compute(concentration, numP1, rMF, seed, partitionSize, Species, pCCEOrder, rDipole, numExtAv, numIntAv, tauMax, points, Omega, Delta, nPulse, LGDriven, DrivenLines, rhoNV, Sparse, meanField), 1:K)
        end
    else
        signal = [zeros(points) for _ in 1:K]
        if ExecuteInCluster
            for k in 1:K
                signal[k] .= Compute(concentration, numP1, rMF, seed, partitionSize, Species, pCCEOrder, rDipole, numExtAv, numIntAv, tauMax, points, Omega, Delta, nPulse, LGDriven, DrivenLines, rhoNV, Sparse, meanField)
            end
        else
            @showprogress for k in 1:K
                signal[k] .= Compute(concentration, numP1, rMF, seed, partitionSize, Species, pCCEOrder, rDipole, numExtAv, numIntAv, tauMax, points, Omega, Delta, nPulse, LGDriven, DrivenLines, rhoNV, Sparse, meanField)
            end
        end
    end

    return signal
end

end