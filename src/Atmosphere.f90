module Atmosphere
    use Constants
    implicit none
    real :: pressure ! [atm] -- pressure on the current atmospheric level
    real :: pSelf ! [atm] -- partial pressure of the considered gaseous species
    real :: pForeign ! [atm] -- foreign pressure ( p - pSelf)
    real :: temperature ! [K] -- temperauture on the current atmospheric level
    real :: density ! [molecule/(cm^2 * km)] -- such density units are needed for having absorption coefficient in km-1
contains
    subroutine fetchAtmosphericParameters()
        ! TODO: read from atmospheric file later
        ! pressure = 1 ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
        ! density = 0.6994E+20 ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
        
        pressure = 0.44126E-02 ! ~ Venus 80 km level
        density = 0.4782E+17 ! ~ Venus 80 km level (see data/Atmospheres/H2O_117.dat)
        ! temperature = refTemperature
        temperature = 203.
        
        pSelf = density * 10. / LOSCHMIDT * temperature/stTemperature
        ! pSelf = 3.8113994E-5
        pForeign = pressure - pSelf
    end subroutine fetchAtmosphericParameters
end module Atmosphere
