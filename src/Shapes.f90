module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy
    implicit none 
contains
    
    ! TODO: add line shape dependency from XX=(nu-nu0) instead of nu
    
    real function simpleLorentz(nu)
        real(kind=DP), intent(in) :: nu ! [cm-1], (gridWV) -- spectral point in which the total contribution from lines is calculated
        real(kind=DP) :: shiftedLineWV ! [cm-1]
        real(kind=DP) :: HWHM

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        HWHM = lorentzHWHM(pressure, temperature, temperatureDependent=.false.)
        simpleLorentz = HWHM / (pi*((nu-shiftedLineWV)**2 + HWHM**2))
    end function simpleLorentz

    real function simpleDoppler(nu)
        real(kind=DP), intent(in) :: nu ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated
        real(kind=DP) :: shiftedLineWV ! [cm-1]
        real(kind=DP) :: HWHM ! cm-1 -- Doppler HWHM

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        HWHM = dopplerHWHM(lineWV, temperature, molarMass, shiftedLineCenter=.false.)
        simpleDoppler = sqrt(log(2.) / (pi * (HWHM**2))) * exp(-(((nu - shiftedLineWV)**2) * log(2.)) / (HWHM**2))
    end function simpleDoppler

end module Shapes
