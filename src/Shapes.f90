module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy
    implicit none
contains
    
    real function simpleLorentz(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, temperature, temperatureDependent=.false.)
        simpleLorentz = HWHM / (pi*(X**2 + HWHM**2))
    end function simpleLorentz

    real function doppler(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !

        real(kind=DP) :: HWHM ! [cm-1] -- Doppler HWHM
  
        HWHM = dopplerHWHM(lineWV, temperature, molarMass)        
        doppler = sqrt(log(2.) / (pi * HWHM**2)) * exp(-(X/HWHM)**2 * log(2.))
    end function doppler

end module Shapes
