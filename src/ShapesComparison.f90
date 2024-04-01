! module ShapesComparison
!     integer :: sampleJ
!     real :: sampleDopplerHWHM, sampleLorentzHWHM
!     real(kind=DP) :: sampleLineWV
!     real(kind=DP), allocatable :: sampleDopplerArray(:,:)
!     real(kind=DP), allocatable :: sampleLorentzArray(:,:)  
    ! character(len=20), parameter :: sampleOutFile = 'calc/sample101.dat'
    ! real(kind=DP), parameter :: delta = 0.01
    ! real(kind=DP), parameter :: sampleStep = 0.00001
    ! integer :: sampleLen, ii
    ! integer, parameter :: sampleRec = 5766

!     implicit none

!     open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
!     read(hitranFileUnit, rec=sampleRec) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
!                             lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
!     sampleLineWV = lineWV

!     sampleDopplerHWHM = dopplerHWHM(pressure, refTemperature, molarMass)
!     sampleLorentzHWHM = simpleLorentzHWHM(pressure)


!     sampleLen = (2 * delta / sampleStep) + 1
!     allocate(sampleDopplerArray(sampleLen,2))
!     allocate(sampleLorentzArray(sampleLen,2))


!     do ii = 1, sampleLen
!         sampleDopplerArray(ii,1) = sampleLineWV - delta + (ii-1) * sampleStep
!         sampleLorentzArray(ii,1) = sampleDopplerArray(ii,1)
!         sampleLorentzArray(ii,2) = simpleLorentz( sampleLorentzArray(ii,1) )
!         sampleDopplerArray(ii,2) = simpleDoppler( sampleDopplerArray(ii,1) )
!     end do

!     open(sampleOutFileUnit, file=sampleOutFile, status='replace', action='write')
!     do ii = 1, sampleLen
!         write(sampleOutFileUnit, '(F20.6, ", ", E20.14, ", ", E20.14)') &
!         sampleDopplerArray(ii,1), sampleLorentzArray(ii,2), &
!         sampleDopplerArray(ii,2)
!     end do
!     close(sampleOutFileUnit)
!     deallocate(sampleDopplerArray)
!     deallocate(sampleLorentzArray)
!     close(hitranFileUnit)
! end module ShapesComparison