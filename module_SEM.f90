!#############################################################################
      module module_SEM
use, intrinsic :: iso_fortran_env, only: dp => real64
!#############################################################################
          real(dp),DIMENSION(:,:,:), ALLOCATABLE :: Vsem ,Usem
          real(dp),DIMENSION(:,:), ALLOCATABLE :: X_EDDY,EPSILO,MOLT
          real(dp),DIMENSION(:,:), ALLOCATABLE :: SIGMA
          real(dp),DIMENSION(:), ALLOCATABLE ::  Ksem
          real(dp),DIMENSION(:) :: X_POINT(3),REYNOLDS(6)
          real(dp),DIMENSION(:) :: TEMP(3),TEMP2(3)
          real(dp), DIMENSION(:,:) ::  R(3,3)
          CHARACTER(len=44) :: FILEGLOBAL
          INTEGER,ALLOCATABLE:: elemyst(:),elemyen(:),elemzst(:),elemzen(:)
          INTEGER,ALLOCATABLE:: iddom(:),ljdom(:),lkdom(:)
      end module module_SEM
