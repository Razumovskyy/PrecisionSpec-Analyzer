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
                                    absorptionCoeffCalc(gridWV, refLineIntensity, shapeFuncPtr)
            end if
            j = j + 1
        spectra(i, 2) = gridAbsorptionCoeff
        end do
    end do

    close(hitranFileUnit)
    
    !! writes spectra array (wavenumber, absorption coeff) to the output ASCII file
    call generateOutput
contains

    real function absorptionCoeffCalc(nu, intensity, lineShape)
        real, intent(in) :: intensity ! [cm-1/(molecule*cm-2)] (refLineIntensity) -- the spectral line intensity for a single molecule per unit volume.
        real(kind=DP), intent(in) :: nu ! [cm-1],  gridWV -- spectral point where absorption coefficitent will be calculated
        procedure(shape), pointer, intent(in) :: lineShape ! [cm] -- normalized line shape function from the MyShapes module

        real(kind=DP) :: X
        real(kind=DP) :: shiftedLineWV

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        X = abs(nu - shiftedLineWV)
        absorptionCoeffCalc = intensity * density * lineShape(X)
    end function absorptionCoeffCalc

end program SpecAnalyzer
