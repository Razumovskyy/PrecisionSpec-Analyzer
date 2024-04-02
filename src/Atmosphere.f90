module Atmosphere
    use Constants
    implicit none
    real :: pressure ! [atm] -- pressure on the current atmospheric level
    real :: temperature ! [K] -- temperauture on the current atmospheric level
    real :: density ! [molecule/(cm^2 * km)] -- such density units are needed for having absorption coefficient in km-1
contains
    subroutine fetchAtmosphericParameters()
        ! read from atmospheric file later
        pressure = 1. ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
        density = 0.6994E+20 ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
        temperature = refTemperature
    end subroutine fetchAtmosphericParameters
end module Atmosphere
