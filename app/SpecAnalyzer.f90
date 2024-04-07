program SpecAnalyzer
    use Constants
    use IO
    use Atmosphere
    use Spectroscopy
    use ShapeFuncInterface
    use Shapes
    implicit none
    
    integer, parameter :: hitranFileUnit = 7777

    real(kind=DP) :: gridWV ! [cm-1] -- more human-readable variable for the x-values of the spectra array
    real :: gridAbsorptionCoeff ! [km-1] -- more human readable variable for the y-values of the spectra array
    
    integer :: startRec ! HITRAN direct access file record number, which corresponds to the first line falling in the [startWV, endWV] interval
    integer :: i, j ! loop variables
    ! ---------------------------------------------------------------------------- !   

    !! reads input from the inputConfig.ini
    call readInputParameters

    !! specifies line shape fucntion to use in calculations
    call fetchLineShapeFunction

    !! reads TIPS from the input file
    call readTIPS
    
    !! fetches pressure, density, temperature and molar mass
    call fetchAtmosphericParameters
    call fetchMolecularParameters

    !! reading from HITRAN direct access file to locate the first spectral line falling in the interval (lineWV, startRec)
    startRec = 1
    do
        open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
        read(hitranFileUnit, rec=startRec) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                        lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
        if (lineWV > startWV - lineCutOff) exit
        startRec = startRec + 1
    end do

    !! allocates resulting 2D array: (wavunumber (cm-1), absorptionCoeff (km-1)) 
    call allocateSpectraArray

    do i = 1, arrayLen
        write(*,*) i, ' of ', arrayLen, ' is processed'
        spectra(i, 1) = startWV + (i-1) * calcPrecision
        gridAbsorptionCoeff = 0.
        gridWV = spectra(i, 1)
        j = startRec
        do
            open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
            read(hitranFileUnit, rec=j) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                    lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign 
            if (lineWV >= endWV + lineCutOff) exit
            if (abs(lineWV - gridWV) < lineCutOff) then
                gridAbsorptionCoeff = gridAbsorptionCoeff + &
                                    absorptionCoeffCalc(gridWV, temperature, shapeFuncPtr)
            end if
            j = j + 1
        spectra(i, 2) = gridAbsorptionCoeff
        end do
    end do

    close(hitranFileUnit)
    
    !! writes spectra array (wavenumber, absorption coeff) to the output ASCII file
    call generateOutput
contains

    real function absorptionCoeffCalc(nu, temperatureParameter, lineShape)
        real, intent(in) :: temperatureParameter ! [K] -- temperature at the current atmospheric level
        real(kind=DP), intent(in) :: nu ! [cm-1],  gridWV -- spectral point where absorption coefficitent will be calculated
        procedure(shape), pointer, intent(in) :: lineShape ! [cm] -- normalized line shape function from the MyShapes module
        ! --------------------------------------------------- !
        real(kind=DP) :: X
        real(kind=DP) :: shiftedLineWV
        real :: intensity ! [cm-1/(molecule*cm-2)] (refLineIntensity) -- the spectral line intensity for a single molecule per unit volume.
        real :: TIPSFactor, boltzmannFactor, emissionFactor
        real :: TIPSOfT
        real :: TIPSOfRefT
        integer :: NTAB_G
        real :: C_G1, C_G2
        real :: t_G1
        integer :: isotopeNum

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        X = abs(nu - shiftedLineWV)

        NTAB_G = (temperatureParameter - 20.0) / 2 + 1
        t_G1 = NTAB_G * 2.0 + 18.
        C_G2 = (temperatureParameter - t_G1)/2.
        C_G1 = 1. - C_G2
        TIPSOfT = C_G1 * TIPS(isotopeNum, NTAB_G) + C_G2 * TIPS(isotopeNum, NTAB_G+1)
        TIPSOfRefT = TIPS(isotopeNum, 139)

        TIPSFactor = TIPSOfT / TIPSOfRefT
        boltzmannFactor = exp(-C2*lineLowerState/temperatureParameter) / exp(-C2*lineLowerState/refTemperature)
        emissionFactor = (1 - exp(-C2*lineWV/temperatureParameter)) / (1 - exp(-C2*lineWV/refTemperature))

        intensity = refLineIntensity * TIPSFactor * boltzmannFactor * emissionFactor
        absorptionCoeffCalc = intensity * density * lineShape(X)
    end function absorptionCoeffCalc

end program SpecAnalyzer
