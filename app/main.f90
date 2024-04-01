program MySpectra
    use MyConstants
    use IO
    use MyAtmosphere
    use MySpectroscopy
    use MyShapeFuncInterface
    use MyShapes
    implicit none
    
    integer, parameter :: hitranFileUnit = 7777
    integer, parameter :: sampleOutFileUnit = 9999

    procedure(shape), pointer :: shapeFuncPtr ! pointer for implementing different line shapes functions

    real(kind=DP) :: gridWV ! [cm-1] -- more human-readable variable for the x-values of the spectra array
    real :: gridAbsorptionCoeff ! [km-1] -- more human readable variable for the y-values of the spectra array
    integer :: startRec ! HITRAN direct access file record number, which corresponds to the first line falling in the [startWV, endWV] interval
    integer :: i, j ! loop variables
    ! ---------------------------------------------------------------------------- !   

    call readInputParams
    call fetchAtmosphericParameters
    call fetchMolecularParameters

    startRec = 1
    do
        open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
        read(hitranFileUnit, rec=startRec) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                        lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
        if (lineWV > startWV - lineCutOff) exit
        startRec = startRec + 1
    end do

    call allocateOutputSpectraArray

    shapeFuncPtr => simpleDoppler
    
    do i = 1, arrayLen
        write(*,*) i, ' of ', arrayLen, ' is processed'
        spectra(i, 1) = startWV + (i-1) * step
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
    
    call generateOutput
contains

    real function absorptionCoeffCalc(nu, intensity, lineShape)
        real, intent(in) :: intensity ! [cm-1/(molecule*cm-2)] (refLineIntensity) -- the spectral line intensity for a single molecule per unit volume.
        real(kind=DP), intent(in) :: nu ! [cm-1],  gridWV -- spectral point where absorption coefficitent will be calculated
        procedure(shape), pointer, intent(in) :: lineShape ! [cm] -- lineShape Function from the MyShapes module

        absorptionCoeffCalc = intensity * density * lineShape(nu)
    end function absorptionCoeffCalc

end program MySpectra


