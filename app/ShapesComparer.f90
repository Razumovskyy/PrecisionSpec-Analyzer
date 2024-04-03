program ShapesComparer
    use Constants
    use CompareIO
    use Atmosphere
    use Spectroscopy
    use ShapeFuncInterface
    use Shapes
    implicit none

    integer, parameter :: hitranFileUnit = 7777
    procedure(shape), pointer :: shape1FuncPtr, shape2FuncPtr

    integer :: loopRecNum
    real(kind=DP) :: baselineWV
    integer :: linesCount
    integer :: baseRecNum
    integer :: i
    ! ---------------------------------------------------------------------------- !  

    !! reads input from the inputConfig.ini
    call readInputParameters
    ! write(*,*) outputShapeFile
    ! pause

    shape1FuncPtr => getLineShapeFunction(lineShape1FuncName)
    shape2FuncPtr => getLineShapeFunction(lineShape2FuncName)

    !! fetches pressure, density, temperature and molar mass
    call fetchAtmosphericParameters
    call fetchMolecularParameters

    !! reading from HITRAN direct access file to locate the first spectral line falling in the interval (lineWV, startRec)
    loopRecNum = 1
    linesCount = 0
    ! write(*,*) lineWV, startWV, endWV
    ! pause
    do
        open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
        read(hitranFileUnit, rec=loopRecNum) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                        lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
        if ((lineWV > startWV) .and. (lineWV < endWV)) then
            baselineWV = lineWV
            baseRecNum = loopRecNum
            linesCount = linesCount + 1
        end if
        if (linesCount > 1) then
            write(*,*) 'More than one lines found in the interval. & 
                        Exited the loop and selected the first line'
            exit
        end if
        if (lineWV > endWV) exit
        loopRecNum = loopRecNum + 1
    end do

    write(*,*) 'baselineWV: ', baselineWV
    write(*,*) 'baserecNum: ', baseRecNum
    ! pause

    call allocateLineShapesArray((baselineWV-lineCutOff), (baselineWV+lineCutOff))

    read(hitranFileUnit, rec=baseRecNum) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                    lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
    close(hitranFileUnit)

    do i = 1, arrayLen
        ! write(*,*) i, ' of ', arrayLen, ' is processed'
        lineShapes(i, 1) = baselineWV - lineCutOff + (i-1) * calcPrecision
        lineShapes(i, 2) = shape1FuncPtr(lineShapes(i,1))
        lineShapes(i, 3) = shape2FuncPtr(lineShapes(i,1))
    end do

    call generateOutput
end program ShapesComparer