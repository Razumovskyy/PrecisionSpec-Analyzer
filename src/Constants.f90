module Constants
    implicit none
    integer, parameter :: DP = selected_real_kind(15,  307) ! selected kind for double precision
    real, parameter :: pi = 3.14159

    ! CGS system is used
    real, parameter :: SPL = 2.99792458e10 ! [cm/s] -- speed of light
    real, parameter :: BOL = 1.3806503e-16 ! [erg/K] -- Boltzmann constant
    real, parameter :: AVOGADRO = 6.02214076e23
    real, parameter :: PLANCK = 6.626070e-27 ! [erg*s]
    real, parameter :: LOSCHMIDT = 2.6867811e25 ! [1/m^3]
    real, parameter :: stTemperature = 273.15 ! [K]
    ! ----------------------------------------------------------------- !
    real, parameter :: C2 = (PLANCK * SPL / BOL) ! [cm*K] -- second radiational constant
    ! ----------------------------------------------------------------- !
    real, parameter :: dopplerCONST = sqrt(2*AVOGADRO*BOL*log(2.)) / SPL
    real, parameter :: refTemperature = 296. ! [K] -- reference temperature for HITRAN data
    real, parameter :: sqln2 = sqrt(log(2.))
end module Constants
