module CompareIO
    use Constants
    use ShapeFuncInterface
    use Shapes
    implicit none
    integer, parameter :: inputConfigUnit = 4443
    character(len=29), parameter :: inputConfigFile = 'shapesConfig.ini'
    integer, parameter :: outputFileUnit = 5558

    ! other variables
    integer :: arrayLen ! length of the array for the line shape functions
    integer :: k ! loop variable

    ! ---------------------------------------------------------------------------- !
    ! INPUT PARAMETERS !
    character(len=20) :: hitranFile
    character(len=30) :: outputShapeFile
    real(kind=DP) :: startWV, endWV ! [cm-1] -- boundaries and grid calcPrecision for the spectral interval for calculating spectra
    real(kind=DP) :: lineCutOff, calcPrecision ! [cm-1] -- line cut-off (measured from the line center)
    character(len=30) :: lineShape1FuncName ! name of the 1st line shape function (custom or standard)
    character(len=30) :: lineShape2FuncName ! name of the 2nd line shape function (custom or standard)
    ! ---------------------------------------------------------------------------- !   
    ! OUTPUT PARAMETERS !
    real(kind=DP), allocatable :: lineShapes(:,:) ! nx3 array: wavenumber (cm-1), lineshape1 (cm), lineshape2 (cm)
    ! ---------------------------------------------------------------------------- !   
contains

    subroutine readInputParameters()
        open(unit=inputConfigUnit, file=inputConfigFile)
        read(inputConfigUnit, '(A20)') hitranFile
        read(inputConfigUnit, '(A30)') outputShapeFile
        read(inputConfigUnit, *) startWV, endWV
        read(inputConfigUnit, *) lineCutOff 
        read(inputConfigUnit, *) calcPrecision
        read(inputConfigUnit, '(A30)') lineShape1FuncName
        read(inputConfigUnit, '(A30)') lineShape2FuncName
    end subroutine readInputParameters

    function getLineShapeFunction(name) result(funcPtr)
        character(len=*), intent(in) :: name
        procedure(shape), pointer :: funcPtr
    
        select case(trim(adjustl(name)))
        case ('simpleLorentz')
            funcPtr => simpleLorentz
        case ('simpleDoppler')
            funcPtr => doppler
        case default
            funcPtr => null()
        end select
    end function getLineShapeFunction

    subroutine allocateLineShapesArray(leftBoundary, rightBoundary)
        real(kind=DP) :: leftBoundary, rightBoundary

        arrayLen = int((rightBoundary-leftBoundary) / calcPrecision) + 1
        allocate(lineShapes(arrayLen, 3))
    end subroutine allocateLineShapesArray

    subroutine generateOutput()
        open(outputFileUnit, file=outputShapeFile)

        do k = 1, arrayLen
            write(outputFileUnit, '(F10.5, ", ", E20.14, ", ", E20.14)') &
            lineShapes(k,1), lineShapes(k,2), lineShapes(k,3)
        end do

        close(outputFileUnit)
        deallocate(lineShapes)
        write(*,*) 'Success !!!'
    end subroutine generateOutput
end module CompareIO