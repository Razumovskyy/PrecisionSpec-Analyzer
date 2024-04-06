module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy
    implicit none
contains
    
    real function simpleLorentz(X)
        !! 
        ! Simplistic Lorentz line shape function in which there is no account for self-broadening and no temperature dependence
        !!
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, includeGammaSelf=.false., includeTemperature=.false.)
        simpleLorentz = HWHM / (pi*(X**2 + HWHM**2))
    end function simpleLorentz

    real function selfSimpleLorentz(X)
        !! 
        ! Lorentz line shape function in which there is no temperature dependence (T = 296 K), but self-broadening is taken into account
        !!
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, includeGammaSelf=.true., includeTemperature=.false.)
        selfSimpleLorentz = HWHM / (pi*(X**2 + HWHM**2))
    end function selfSimpleLorentz

    real function noSelfLorentz(X)
        !! 
        ! Temperature-dependent Lorentz line shape with no account for self-broadening
        !!
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, includeGammaSelf=.true., includeTemperature=.true., temperatureParameter=temperature)
        noSelfLorentz = HWHM / (pi*(X**2 + HWHM**2))
    end function noSelfLorentz

    real function lorentz(X)
        !! 
        ! Temperature-dependent Lorentz line shape with no account for self-broadening
        !!
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, includeGammaSelf=.true., includeTemperature=.true., temperatureParameter=temperature)
        lorentz = HWHM / (pi*(X**2 + HWHM**2))
    end function lorentz

    real function doppler(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM ! [cm-1] -- Doppler HWHM
  
        HWHM = dopplerHWHM(lineWV, temperature, molarMass)        
        doppler = sqrt(log(2.) / (pi * HWHM**2)) * exp(-(X/HWHM)**2 * log(2.))
    end function doppler

end module Shapes
