PROGRAM accumulated_stats

!=======================================================================
! Purpose:
! --------
! *accumulated_stats* Read an accumulated field from GRIB output files 
! from a whole year, run and write statistical data for these.
!
! Interface:
! ----------
! Variable     Type         Content
! ========     ====         =======
!
! Author:      Kristian P. Nielsen,                           2021-06-07
! -------
!
! Modifications:
! --------------
!
!=======================================================================
! 1.    Declarations:
!=======================================================================

USE ECCODES

IMPLICIT NONE

INTEGER, PARAMETER              :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER              :: JPRM = SELECTED_REAL_KIND(6,37)

!-----------------------------------------------------------------------
! 1.1   LOCAL VARIABLES:
!-----------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER   :: harmkx=1069
INTEGER(KIND=JPIM), PARAMETER   :: harmky=1269
REAL(KIND=JPRM),    PARAMETER   :: pixel_area=6250000 ! m²
REAL(KIND=JPRM),    PARAMETER   :: trillion=1.E12 
REAL(KIND=JPRM),    PARAMETER   :: zero=0.0
CHARACTER(LEN=2),   PARAMETER   :: shortName1='tp'
CHARACTER(LEN=2),   PARAMETER   :: shortName2='mn2t'
CHARACTER(LEN=2),   PARAMETER   :: shortName3='mx2t'
CHARACTER(LEN=47)               :: dayfile
CHARACTER(LEN=16)               :: fstem
CHARACTER(LEN=49)               :: harmfile
CHARACTER(LEN=5)                :: hour_int
CHARACTER(LEN=15)               :: harmdir
CHARACTER(LEN=27)               :: igjfile
CHARACTER(LEN=2)                :: ini_hour
CHARACTER(LEN=13)               :: lonlats
CHARACTER(LEN=60)               :: sname
CHARACTER(LEN=3)                :: strday
CHARACTER(LEN=4)                :: stryear
CHARACTER(LEN=46)               :: yearfile
INTEGER(KIND=JPIM)              :: day
INTEGER(KIND=JPIM)              :: h
INTEGER(KIND=JPIM)              :: i
INTEGER(KIND=JPIM)              :: j
INTEGER(KIND=JPIM)              :: jcar
INTEGER(KIND=JPIM)              :: jmir
INTEGER(KIND=JPIM)              :: ierr
INTEGER(KIND=JPIM)              :: igrib
INTEGER(KIND=JPIM)              :: iunit
INTEGER(KIND=JPIM)              :: nday
INTEGER(KIND=JPIM)              :: numberofpoints
INTEGER(KIND=JPIM)              :: ogrib
INTEGER(KIND=JPIM)              :: ounit
INTEGER(KIND=JPIM)              :: t
INTEGER(KIND=JPIM)              :: year
REAL(KIND=JPRM)                 :: daily_prec(harmkx*harmky,366)
REAL(KIND=JPRM)                 :: GRL_glacier_mask(harmkx*harmky)
REAL(KIND=JPRM)                 :: lat
REAL(KIND=JPRM)                 :: long
REAL(KIND=JPRM)                 :: ISL_glacier_mask(harmkx*harmky)
REAL(KIND=JPRM)                 :: JM_glacier_mask(harmkx*harmky)
REAL(KIND=JPRM)                 :: tmpmask
REAL(KIND=JPRM)                 :: yearly_GRL_glacier_prec
REAL(KIND=JPRM)                 :: yearly_ISL_glacier_prec
REAL(KIND=JPRM)                 :: yearly_JM_glacier_prec
REAL(KIND=JPRM)                 :: yearly_prec(harmkx*harmky)
REAL(KIND=JPRM),DIMENSION(:),ALLOCATABLE :: values

!=======================================================================
! 2.    Basic definitions:
!=======================================================================

READ(5,*)year ! Streaming input $YYYY
harmdir = '/data/kpn/CARRA'
fstem   = 'CARRA-West_prec_'
WRITE(stryear,'(I4)')year
nday    = 365
IF (MOD(year,4).EQ.0)nday=366

!-----------------------------------------------------------------------
! 2.1   Reading the glacier masks:
!-----------------------------------------------------------------------

GRL_glacier_mask = zero
ISL_glacier_mask = zero
JM_glacier_mask  = zero
lonlats = 'CW_coords.dat'
igjfile = 'glacier_mask_GRL_ISL_JM.txt'
OPEN(1,FILE=harmdir//'/'//lonlats)
OPEN(2,FILE=harmdir//'/'//igjfile)
DO j=1,harmky
  DO i=1,harmkx
    READ(1,102)long,lat
    READ(2,101)tmpmask
    IF (lat.LT.68.0.AND.long.GT.-27.0) THEN
      ISL_glacier_mask(i+harmkx*(j-1)) = tmpmask
    ELSE IF(lat.GT.70.0.AND.lat.LT.75.0.AND.long.GT.-15.0) THEN
      JM_glacier_mask(i+harmkx*(j-1))  = tmpmask
    ELSE
      GRL_glacier_mask(i+harmkx*(j-1)) = tmpmask
    ENDIF
  ENDDO
  READ(1,*)
  READ(2,*)
ENDDO
CLOSE(1)
CLOSE(2)
!WRITE(*,*)'GRL glacier area [m²]: ',SUM(GRL_glacier_mask)*pixel_area
!WRITE(*,*)'ISL glacier area [m²]: ',SUM(ISL_glacier_mask)*pixel_area
!WRITE(*,*)'JM glacier area [m²]: ',SUM(JM_glacier_mask)*pixel_area

!=======================================================================
! 3.    Cycle through, read and sum accumulated data for 06-06 UTC runs:
!=======================================================================

daily_prec              = zero
yearly_prec             = zero
yearly_GRL_glacier_prec = zero
yearly_ISL_glacier_prec = zero
yearly_JM_glacier_prec  = zero
DO t=0,12,12
  WRITE(ini_hour,'(I2.2)')t
  DO h=1,4
    IF ( h.eq.1 ) hour_int='9-6  '
    IF ( h.eq.2 ) hour_int='12-9 '
    IF ( h.eq.3 ) hour_int='15-12'
    IF ( h.eq.4 ) hour_int='18-15'
    harmfile = fstem//stryear//'_3h_acc_'//ini_hour//'_UTC_fl_'// &
     &         TRIM(hour_int)//'h.grib'
    CALL codes_open_file(iunit,harmdir//'/'//harmfile,'r')
    DO day=1,nday
      WRITE(strday,'(I3.3)')day
      CALL codes_grib_new_from_file(iunit,igrib,ierr)
      IF (day.EQ.1)CALL codes_clone(igrib,ogrib)
      CALL codes_get(igrib,'numberOfPoints',numberofpoints)
      IF (harmkx*harmky.NE.numberofpoints) THEN
        WRITE(*,*)'HARMONIE data array size error!'
      ENDIF
      CALL codes_get(igrib,'shortName',sname)
      IF (sname(1:3).EQ.'tp ') THEN
        ALLOCATE(values(numberofpoints))
        CALL codes_get(igrib,"values",values)
        daily_prec(:,day) = daily_prec(:,day) + values
        yearly_prec       = yearly_prec + values
        yearly_GRL_glacier_prec = yearly_GRL_glacier_prec + &
         &                        SUM(values*GRL_glacier_mask)* &
         &                        pixel_area
        yearly_ISL_glacier_prec = yearly_ISL_glacier_prec + &
         &                        SUM(values*ISL_glacier_mask)* &
         &                        pixel_area
        yearly_JM_glacier_prec  = yearly_JM_glacier_prec + &
         &                        SUM(values*JM_glacier_mask)* &
         &                        pixel_area
        DEALLOCATE(values)
      ENDIF
      CALL codes_release(igrib)
    ENDDO
    CALL codes_close_file(iunit)
  ENDDO
ENDDO
yearly_GRL_glacier_prec = yearly_GRL_glacier_prec / trillion ! -> Gt
yearly_ISL_glacier_prec = yearly_ISL_glacier_prec / trillion ! -> Gt
yearly_JM_glacier_prec  = yearly_JM_glacier_prec / trillion ! -> Gt
WRITE(*,204)'GRL glacier ',stryear,' 06-06 UTC total precipitation [Gt]: ', &
 &          yearly_GRL_glacier_prec
WRITE(*,304)'ISL glacier ',stryear,' 06-06 UTC total precipitation [Gt]: ', &
 &          yearly_ISL_glacier_prec
WRITE(*,404)'JM glacier ',stryear,' 06-06 UTC total precipitation [Gt]: ', &
 &          yearly_JM_glacier_prec

!!-----------------------------------------------------------------------
!! 3.1   Write daily accumulated precipitation data:
!!-----------------------------------------------------------------------
!
!DO day=1,nday
  !WRITE(strday,'(I3.3)')day
  !CALL codes_set(ogrib,'level',0)
  !CALL codes_set(ogrib,'stepType','accum')
  !CALL codes_set(ogrib,'stepUnits','D')   ! Daily accumulated values
  !CALL codes_set(ogrib,'stepRange','0-1')
  !CALL codes_set(ogrib,'shortName',shortName1)
  !ALLOCATE(values(numberofpoints))
  !values = daily_prec(:,day)
  !CALL codes_set(ogrib,'values',values)
  !DEALLOCATE(values)
  !dayfile = fstem//stryear//'_'//strday//'_24h_acc_06-06_UTC.grib'
  !CALL codes_open_file(ounit,harmdir//'/'//dayfile,'w')
  !CALL codes_write(ogrib,ounit)
  !CALL codes_close_file(ounit)
!ENDDO

!-----------------------------------------------------------------------
! 3.2   Write yearly accumulated precipitation data:
!-----------------------------------------------------------------------

CALL codes_set(ogrib,'level',0)
CALL codes_set(ogrib,'stepType','accum')
CALL codes_set(ogrib,'stepUnits','Y')   ! Yearly accumulated values
CALL codes_set(ogrib,'stepRange','0-1')
CALL codes_set(ogrib,'shortName',shortName1)
ALLOCATE(values(numberofpoints))
values = yearly_prec
CALL codes_set(ogrib,'values',values)
DEALLOCATE(values)
yearfile = fstem//stryear//'_yearly_acc_06-06_UTC.grib'
CALL codes_open_file(ounit,harmdir//'/'//yearfile,'w')
CALL codes_write(ogrib,ounit)
CALL codes_close_file(ounit)

!=======================================================================
! 4.    Cycle through, read and sum accumulated data for 03-03 UTC runs:
!=======================================================================

daily_prec              = zero
yearly_prec             = zero
yearly_GRL_glacier_prec = zero
yearly_ISL_glacier_prec = zero
yearly_JM_glacier_prec  = zero
DO t=0,12,12
  WRITE(ini_hour,'(I2.2)')t
  DO h=1,4
    IF ( h.eq.1 ) hour_int='6-3  '
    IF ( h.eq.2 ) hour_int='9-6  '
    IF ( h.eq.3 ) hour_int='12-9 '
    IF ( h.eq.4 ) hour_int='15-12'
    harmfile = fstem//stryear//'_3h_acc_'//ini_hour//'_UTC_fl_'// &
     &         TRIM(hour_int)//'h.grib'
    CALL codes_open_file(iunit,harmdir//'/'//harmfile,'r')
    DO day=1,nday
      WRITE(strday,'(I3.3)')day
      CALL codes_grib_new_from_file(iunit,igrib,ierr)
      IF (day.EQ.1)CALL codes_clone(igrib,ogrib)
      CALL codes_get(igrib,'numberOfPoints',numberofpoints)
      IF (harmkx*harmky.NE.numberofpoints) THEN
        WRITE(*,*)'HARMONIE data array size error!'
      ENDIF
      CALL codes_get(igrib,'shortName',sname)
      IF (sname(1:3).EQ.'tp ') THEN
        ALLOCATE(values(numberofpoints))
        CALL codes_get(igrib,"values",values)
        daily_prec(:,day) = daily_prec(:,day) + values
        yearly_prec       = yearly_prec + values
        yearly_GRL_glacier_prec = yearly_GRL_glacier_prec + &
         &                        SUM(values*GRL_glacier_mask)* &
         &                        pixel_area
        yearly_ISL_glacier_prec = yearly_ISL_glacier_prec + &
         &                        SUM(values*ISL_glacier_mask)* &
         &                        pixel_area
        yearly_JM_glacier_prec  = yearly_JM_glacier_prec + &
         &                        SUM(values*JM_glacier_mask)* &
         &                        pixel_area
        DEALLOCATE(values)
      ENDIF
      CALL codes_release(igrib)
    ENDDO
    CALL codes_close_file(iunit)
  ENDDO
ENDDO

!-----------------------------------------------------------------------
! 4.1   Write yearly accumulated precipitation data:
!-----------------------------------------------------------------------

CALL codes_set(ogrib,'level',0)
CALL codes_set(ogrib,'stepType','accum')
CALL codes_set(ogrib,'stepUnits','Y') ! Yearly accumulated values
CALL codes_set(ogrib,'stepRange','0-1')
CALL codes_set(ogrib,'shortName',shortName1)
ALLOCATE(values(numberofpoints))
values = yearly_prec
CALL codes_set(ogrib,'values',values)
DEALLOCATE(values)
yearfile = fstem//stryear//'_yearly_acc_03-03_UTC.grib'
CALL codes_open_file(ounit,harmdir//'/'//yearfile,'w')
CALL codes_write(ogrib,ounit)
CALL codes_close_file(ounit)
yearly_GRL_glacier_prec = yearly_GRL_glacier_prec / trillion ! -> Gt
yearly_ISL_glacier_prec = yearly_ISL_glacier_prec / trillion ! -> Gt
yearly_JM_glacier_prec  = yearly_JM_glacier_prec / trillion ! -> Gt
WRITE(*,204)'GRL glacier ',stryear,' 03-03 UTC total precipitation [Gt]: ', &
 &          yearly_GRL_glacier_prec
WRITE(*,304)'ISL glacier ',stryear,' 03-03 UTC total precipitation [Gt]: ', &
 &          yearly_ISL_glacier_prec
WRITE(*,404)'JM glacier ',stryear,' 03-03 UTC total precipitation [Gt]: ', &
 &          yearly_JM_glacier_prec

101 FORMAT(F6.2)
102 FORMAT(2F8.3)
204 FORMAT(A12,A4,A37,F6.1)
304 FORMAT(A12,A4,A37,F7.2)
404 FORMAT(A12,A4,A37,F8.3)

END
