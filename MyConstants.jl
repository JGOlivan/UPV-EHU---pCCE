#############################################################################################################
################                    This module contains necessary constants                 ################
#############################################################################################################

module MyConstants

    export hbar, mu0, a, ge, sigmax, sigmay, sigmaz, sx, sy, sz

    # Everyday constants

    const hbar      = 1.054571817e-34       # Dirac constant (J·s)         
    const mu0       = 1.256637061e-6        # Magnetic permeability of vacuum (N/A^2)             
    const a         = 3.567e-10             # Diamond Lattice constant (m)       
    const ge        = -28.0249513861e9      # Electron Gyromagnetic Ratio (Hz/T)    

    # ------------
    # OPERATORS
    # ------------

    # Pauli operators

    const sigmax = [0 1; 1 0]
    const sigmay = [0 -im; im 0]
    const sigmaz = [1 0; 0 -1]

    # Spin 1 operators

    const sx = 1/sqrt(2) * [0 1 0; 1 0 1; 0 1 0]
    const sy = 1/sqrt(2) * [0 -im 0; im 0 -im; 0 im 0]
    const sz = [1 0 0; 0 0 0; 0 0 -1]

end