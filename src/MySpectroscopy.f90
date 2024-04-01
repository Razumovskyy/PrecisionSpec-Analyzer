module MySpectroscopy
    use MyConstants, only: DP, pi, dopplerCONST
    use MyAtmosphere
    implicit none
    real(kind=DP) :: lineWV ! [cm-1] -- current spectral line
    real :: refLineIntensity ! [cm-1/(molecule*cm-2)] -- spectral line intensity at refTemperature=296 K
    real :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM at refTemperature=296 K
    real :: gammaSelf ! [cm-1/atm] -- self-broadened component of Lorentz HWHM
    real :: lineLowerState ! [cm-1] -- lower state energy of the transition
    real :: foreignTempCoeff ! [dimensionless] (coefficient for temperature dependence of gammaForeign)
    integer :: jointMolIso ! [dimensionless] (joined reference to Molecule number (MOL) and Isotopologue number (ISO))
    real :: deltaForeign ! [cm-1/atm] (pressure shift of the line position at 296 K and 1 atm)
    real :: molarMass ! [g/mol] -- current species molar mass

contains
 
    subroutine fetchMolecularParameters()
        molarMass = 18. 
    end subroutine fetchMolecularParameters
    
    real(kind=DP) function shiftedLinePosition(pressureParameter)
        real :: pressureParameter

        shiftedLinePosition = lineWV + deltaForeign * pressureParameter
    end function shiftedLinePosition

    real function dopplerHWHM(pressureParameter, temperatureParameter, molarMassParameter)
        real, intent(in) :: pressureParameter
        real, intent(in) :: temperatureParameter
        real, intent(in) :: molarMassParameter

        dopplerHWHM = dopplerCONST * shiftedLinePosition(pressureParameter) * sqrt(temperatureParameter / molarMassParameter) 
    end function DopplerHWHM

    real function simpleLorentzHWHM(pressureParameter)
        real, intent(in) :: pressureParameter

        simpleLorentzHWHM = gammaForeign * pressureParameter
    end function simpleLorentzHWHM

end module MySpectroscopy
