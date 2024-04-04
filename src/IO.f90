module IO
    use Constants
    use ShapeFuncInterface
    use Shapes    
    implicit none
    integer, parameter :: inputConfigUnit = 6663
    character(len=29), parameter :: inputConfigFile = 'inputConfig.ini'
    
    integer, parameter :: outputFileUnit = 7778
    
    procedure(shape), pointer :: shapeFuncPtr ! pointer for implementing different line shapes functions

    integer :: arrayLen ! length of the spectra array
    integer :: k ! loop variable

    ! ---------------------------------------------------------------------------- !
    ! INPUT PARAMETERS !
    character(len=20) :: hitranFile
    character(len=22) :: outputFile
    real(kind=DP) :: startWV, endWV, calcPrecision ! [cm-1] -- boundaries and grid calcPrecision for the spectral interval for calculating spectra 
    real(kind=DP) :: lineCutOff ! [cm-1] -- line cut-off (measured from the line center)
    character(len=30) :: lineShapeFuncName ! name of the line shape function (custom or standard)
    ! ---------------------------------------------------------------------------- !   
    ! OUTPUT PARAMETERS !
    real(kind=DP), allocatable :: spectra(:,:) ! 2D array: first item is grid wavenumber in [cm-1] and the second is the absorption coefficient in [km-1]
    ! ---------------------------------------------------------------------------- !   
contains
    
    ! TODO: add molType to input and automatic fetching needed HITRAN file
    subroutine readInputParameters()
        open(unit=inputConfigUnit, file=inputConfigFile, status='old', action='read')
        read(inputConfigUnit, '(A20)') hitranFile
        read(inputConfigUnit, '(A22)') outputFile
        read(inputConfigUnit, *) startWV, endWV, calcPrecision
        read(inputConfigUnit, *) lineCutOff
        read(inputConfigUnit, '(A30)') lineShapeFuncName
        close(inputConfigUnit)
    end subroutine readInputParameters

    !! If you introduce custom line shape function, you need to add
    !! another case statement in this section
    subroutine fetchLineShapeFunction()
        select case(trim(adjustl(lineShapeFuncName)))
        case ('simpleLorentz')
            shapeFuncPtr => simpleLorentz
        case ('doppler')
            shapeFuncPtr => doppler
        end select
    end subroutine fetchLineShapeFunction

    subroutine allocateSpectraArray()
        arrayLen = int((endWV-startWV) / calcPrecision) + 1
        allocate(spectra(arrayLen, 2))
    end subroutine allocateSpectraArray
    
    subroutine generateOutput()
        open(outputFileUnit, file=outputFile)

        do k = 1, arrayLen
            write(outputFileUnit, '(F8.3, ", ", E20.14)') spectra(k,1), spectra(k, 2)
        end do

        close(outputFileUnit)
        deallocate(spectra)
    end subroutine generateOutput
    
end module IO
