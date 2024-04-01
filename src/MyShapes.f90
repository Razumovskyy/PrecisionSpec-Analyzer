module MyShapes
    use MyConstants
    use MyAtmosphere
    use MySpectroscopy
    implicit none 
contains
    
    real function simpleLorentz(nu)
        real(kind=DP), intent(in) :: nu ! [cm-1], (gridWV) -- spectral point in which the total contribution from lines is calculated
        real(kind=DP) :: shiftedLineWV ! [cm-1]

        !shiftedLineWV = shiftedLinePosition(pressure)
        shiftedLineWV = lineWV
        simpleLorentz = (gammaForeign*pressure) / (pi*((nu-shiftedLineWV)**2 + (gammaForeign*pressure)**2))
    end function simpleLorentz

    real function simpleDoppler(nu)
        real(kind=DP), intent(in) :: nu ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated
        real :: refAlphaD ! cm-1 -- Doppler HWHM at reference temperature of 296 K (refTemperature)
        real(kind=DP) :: shiftedLineWV ! [cm-1]

        shiftedLineWV = shiftedLinePosition(pressure)
        refAlphaD = DopplerHWHM(pressure, refTemperature, molarMass)
        simpleDoppler = sqrt(log(2.) / (pi * (refAlphaD**2))) * exp(-(((nu - shiftedLineWV)**2) * log(2.)) / (refAlphaD**2))
    end function simpleDoppler

end module MyShapes