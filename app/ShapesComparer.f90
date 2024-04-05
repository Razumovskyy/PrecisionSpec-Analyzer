program ShapesComparer
    use Constants
    use CompareIO
    use Atmosphere
    use Spectroscopy
    use ShapeFuncInterface
    use Shapes
    implicit none

    integer, parameter :: hitranFileUnit = 7777
    
    procedure(shape), pointer :: shape1FuncPtr, shape2FuncPtr ! pointers to the line shape functions to be compared
    real(kind=DP) :: baselineWV ! (cm-1) -- center of the considered spectroscopic line
    
    ! ---------------------------------------------------------------------------- !
    character(len=20) :: baseLineWVStr, startWVStr, endWVStr ! string representation of the baseline 
    integer :: linesCount ! number of the liens in the interval [startWV, endWV]
    integer :: i, loopRecNum ! loop variables
    ! ---------------------------------------------------------------------------- !  

    !! reads input from the inputConfig.ini
    call readInputParameters

    !! specifies line shapes functions to be compared
    shape1FuncPtr => getLineShapeFunction(lineShape1FuncName)
    shape2FuncPtr => getLineShapeFunction(lineShape2FuncName)

    !! fetches pressure, density, temperature and molar mass
    call fetchAtmosphericParameters
    call fetchMolecularParameters

    !! reading from HITRAN direct access file to locate the first spectral line falling in the interval (lineWV, startRec)
    loopRecNum = 1
    linesCount = 0

    do
        open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
        read(hitranFileUnit, rec=loopRecNum) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                        lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
        if ((lineWV > startWV) .and. (lineWV < endWV)) then
            linesCount = linesCount + 1
            if (linesCount > 1) then
                write(baseLineWVStr, '(F10.5)') baselineWV
                write(*,'(3A)') 'More than one line found in the interval from the input. ', &
                                'Selected the first line: ' // trim(adjustl(baseLineWVStr)) &
                                // ' cm-1'
                exit
            end if
            baselineWV = lineWV
        end if
        if (lineWV > endWV) exit
        loopRecNum = loopRecNum + 1
    end do

    close(hitranFileUnit)

    if (linesCount == 0) then
        write(startWVStr, '(F10.5)') startWV
        write(endWVStr, '(F10.5)') endWV
        write(*,'(2A)') 'Error: no lines found in the interval: ' // trim(adjustl(startWVStr)), &
                    ' and ' // trim(adjustl(endWVStr))
        error stop
    end if

    call allocateLineShapesArray((baselineWV-lineCutOff), (baselineWV+lineCutOff))

    ! TODO: check that line shape is not calculated twice (left and right wings are the same) -- may be redundant
    do i = 1, arrayLen
        lineShapes(i, 1) = baselineWV - lineCutOff + (i-1) * calcPrecision
        lineShapes(i, 2) = shape1FuncPtr(abs(lineShapes(i,1) - lineWV))
        lineShapes(i, 3) = shape2FuncPtr(abs(lineShapes(i,1) - lineWV))
    end do
    call generateOutput
end program ShapesComparer