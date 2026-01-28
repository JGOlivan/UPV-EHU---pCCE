#############################################################################################################
################                This modules define some useful custom types                 ################
#############################################################################################################

module MyTypes

    using StaticArrays

    using ..MyConstants

    export Position, spinOperators, Spin, Driving

    struct Position{T <: Real} <: FieldVector{3, T}
        x :: T
        y :: T
        z :: T
    end

    struct spinOperators{T <: AbstractArray}

        Sx :: T 
        Sy :: T                 # Central spin operators.
        Sz :: T

        Jx :: Vector{T}
        Jy :: Vector{T}         # Bath spin operators.
        Jz :: Vector{T}

    end

    mutable struct Spin{S <: AbstractString, P <: Position{<: Real}}

        Pos         :: P            # A Position type describing the position of the spin in the lattice.
        ID          :: Int64        # A unique label for the spin.
        NucSpin     :: Float64      # Value of the nuclear spin.
        Species     :: S            # Chemical species of the spin.
        Partition   :: Int64        # For pCCE, a label to the partition to which the spin is assigned.
        Branch      :: Int64        # A reference to the resonant frequency of the spin.
        GyroRatio   :: Float64      # The gyromagnetic ratio of the spin.
        Driven      :: Bool         # True if the spin is being driven, false otherwise.
        LGDriven    :: Bool         # True if the spin is LG driven, false otherwise.

    end

    mutable struct Driving{F <: Real, S <: AbstractString} # Single bath driving tone

        Omega           :: F        # The Rabi strengh of the driving.
        Delta           :: F        # The detuning of the driving.
        alpha           :: F        # The angle of the driving. Has no effect for standard LG.
        LG              :: S        # The type of LG driving (LG, FSLG, LG4, ...). Only LG for now.

    end

end








