module IO
    use MyConstants
    implicit none
    integer, parameter :: inputConfigUnit = 6663
    character(len=29), parameter :: inputConfigFile = 'inputConfig.ini'
    
    integer, parameter :: outputFileUnit = 7778

    integer :: arrayLen ! length of the spectra array
    integer :: k ! loop variable

    ! ---------------------------------------------------------------------------- !
    ! INPUT PARAMETERS !
    character(len=20) :: hitranFile
    character(len=22) :: outputFile
    real(kind=DP) :: startWV, endWV, step ! [cm-1] -- boundaries and grid step for the spectral interval for calculating spectra 
    real(kind=DP) :: lineCutOff ! [cm-1] -- line cut-off (measured from the line center)  
    ! ---------------------------------------------------------------------------- !   
    ! OUTPUT PARAMETERS !
    real(kind=DP), allocatable :: spectra(:,:) ! 2D array: first item is grid wavenumber in [cm-1] and the second is the absorption coefficient in [km-1]
    ! ---------------------------------------------------------------------------- !   
contains
    subroutine readInputParams()
        open(unit=inputConfigUnit, file=inputConfigFile, status='old', action='read')
        read(inputConfigUnit, '(A20)') hitranFile
        read(inputConfigUnit, '(A22)') outputFile
        read(inputConfigUnit, *) startWV, endWV, step
        read(inputConfigUnit, *) lineCutOff
        close(inputConfigUnit)
    end subroutine readInputParams

    subroutine allocateOutputSpectraArray()
        arrayLen = int((endWV-startWV) / step) + 1
        allocate(spectra(arrayLen, 2))
    end subroutine allocateOutputSpectraArray
    
    subroutine generateOutput()
        open(outputFileUnit, file=outputFile, status='replace', action='write')

        do k = 1, arrayLen
            write(outputFileUnit, '(F8.3, ", ", E20.14)') spectra(k,1), spectra(k, 2)
        end do

        close(outputFileUnit)
        deallocate(spectra)
    end subroutine generateOutput
    
end module IO
