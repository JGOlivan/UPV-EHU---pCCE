#############################################################################################################
################                This modules define some useful custom types                 ################
#############################################################################################################

module MyTypes

    using StaticArrays

    using ..MyConstants

    export Position, spinOperators, Spin, Driving

    """
        Position{T <: Real} <: FieldVector{3, T}
    
        A struct with fields x, y, z that represent the cartesian coordinates of a spin. 
    """
    struct Position{T <: Real} <: FieldVector{3, T}
        x :: T
        y :: T
        z :: T
    end

    """
        spinOperators{T <: AbstractArray}

        A struct with fields Sx, Sy, Sz (central spin operators), Jx, Jy, Jz (bath spin operators), Ident (Identity matrix) and sz (z spin operator for a central spin in a reduced subspace; for the particular case of the NV).
    """
    struct spinOperators{T <: AbstractArray, D <: AbstractArray}

        Sx      :: T 
        Sy      :: T                 # Central spin operators.
        Sz      :: T

        Jx      :: Vector{T}
        Jy      :: Vector{T}         # Bath spin operators.
        Jz      :: Vector{T}

        Ident   :: D                 # Identity operator
        sz      :: T                 # NV subspace z-operator

        Upi      :: T                # pi-pulse operator

    end

    """
        Spin{S <: AbstractString, P <: Position{<: Real}}

        A mutable struct that describes a single bath spin. The fields are:
            - Pos           :: A Position type describing the position of the spin in the lattice.
            - ID            :: A unique label for the spin.
            - NucSpin       :: Value of the nuclear spin.
            - Species       :: The chemical species of the spin.
            - Partition     :: For pCCE, a label to the partition to which the spin is asigned.
            - Branch        :: A reference to the resonant frequency of the spin (i.e., the corresponding resonance line).
            - GyroRatio     :: The gyromagnetic ration of the spin.
            - Driven        :: True if the spin is being driven, false otherwise.
            - LGDriven      :: True if the spin is being LG driven, false otherwise. (Can be false when Driven is true but cannot be true is Driven is false).
            - Ax            :: The x-component of the interaction vector with the central spin.
            - Ay            :: The y-component of the interaction vector with the central spin.
            - Az            :: The z-component of the interaction vector with the central spin.
    """
    mutable struct Spin{S <: AbstractString, P <: Position{<: Real}, F <: Real}

        Pos         :: P            
        ID          :: Int64        
        NucSpin     :: F      
        Species     :: S            
        Partition   :: Int64        
        Branch      :: Int64        
        GyroRatio   :: F      
        Driven      :: Bool         
        LGDriven    :: Bool  
        Ax          :: F
        Ay          :: F  
        Az          :: F

    end

    """
        Driving{F <: Real, S <: AbstractString}

        A mutable struct describing a single bath driving tone. The fields are:
            - Omega         :: The Rabi frequency of the driving, i.e., the power. 
            - Delta         :: The detuning of the driving.
            - alpha         :: The angle the driving. For now it does not change a thing as we only consider standard LG.
            - LG            :: The type of LG driving (LG, FSLG, LG4, ...). For now it only accepts LG.
    """
    mutable struct Driving{F <: Real, S <: AbstractString} 

        Omega           :: F        
        Delta           :: F       
        alpha           :: F       
        LG              :: S      

    end

end








