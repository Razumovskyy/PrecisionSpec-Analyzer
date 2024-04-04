module Spectroscopy
    use Constants, only: DP, pi, dopplerCONST, refTemperature
    use Atmosphere
    implicit none
    real(kind=DP) :: lineWV ! [cm-1] -- current spectral line, transition wavenumber
    real :: refLineIntensity ! [cm-1/(molecule*cm-2)] -- spectral line intensity at refTemperature=296 K
    real :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM at refTemperature=296 K
    real :: gammaSelf ! [cm-1/atm] -- self-broadened component of Lorentz HWHM
    real :: lineLowerState ! [cm-1] -- lower state energy of the transition
    real :: foreignTempCoeff ! [dimensionless] (coefficient for temperature dependence of gammaForeign)
    integer :: jointMolIso ! [dimensionless] custom variable: joined reference to Molecule number (MOL) and Isotopologue number (ISO)
    real :: deltaForeign ! [cm-1/atm] (pressure shift of the line position at 296 K and 1 atm)
    real :: molarMass ! [g/mol] -- current species molar mass

contains
 
    subroutine fetchMolecularParameters()
        molarMass = 18. 
    end subroutine fetchMolecularParameters
    
    real(kind=DP) function shiftedLinePosition(lineWVParameter, pressureParameter)
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: pressureParameter

        shiftedLinePosition = lineWVParameter + deltaForeign * pressureParameter
    end function shiftedLinePosition

    real function dopplerHWHM(lineWVParameter, temperatureParameter, molarMassParameter)
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: temperatureParameter
        real, intent(in) :: molarMassParameter

        dopplerHWHM = dopplerCONST * shiftedLinePosition(lineWVParameter, pressure) * &
                sqrt(temperatureParameter / molarMassParameter)

    end function dopplerHWHM

    real function lorentzHWHM(pressureParameter, temperatureParameter, temperatureDependent)
        ! TODO: add partial pressure logic (p_self, p_air)
        real, intent(in) :: pressureParameter
        real, intent(in) :: temperatureParameter
        logical, optional, intent(in) :: temperatureDependent

        logical :: isTemperatureDependent
        isTemperatureDependent = .false.
        if (present(temperatureDependent)) isTemperatureDependent = temperatureDependent

        if (.not. isTemperatureDependent) then
            ! temperature is set to 296 K and partial pressure is not counted
            lorentzHWHM = gammaForeign * pressureParameter
        else
            lorentzHWHM = ((temperatureParameter / refTemperature)**foreignTempCoeff) * (gammaForeign * pressureParameter)
        end if
    end function lorentzHWHM

end module Spectroscopy
