! **********************************************************************
!  Yokota's Undeveloped Model
!  - Sphere Coordinate System
!  - Reduced Gravity
!  - Equatorial Region
!  - Orlanski Radiation condition at Northern and Southern boundaries
!
!           Target Region: 30S-30N, 120E-70W
!           Computational Region : 50S-50N, 100E-60W
!
!                              5th release(orc) 2017/10/10, by h.yokota
!                              4th release(dp) 2017/07/07, by h.yokota
!                              3rd release(f90) 2017/05/31, by h.yokota
!                              2nd release 2003/02/10 by h.yokota
!                              1st release 1995/11/09 by h.yokota
! **********************************************************************

MODULE commons
  IMPLICIT NONE
  INTEGER, PARAMETER :: ix=203, jy=103
  DOUBLE PRECISION :: h(ix,jy,-2:1),u(ix,jy,-2:1),v(ix,jy,-2:1)
  DOUBLE PRECISION :: bc,dt1,dt2,dx,dy,x0,y0
  DOUBLE PRECISION :: ah(jy),ahmn,ahmx,rg,hini,omega,pi,re,rho1,rho2
  DOUBLE PRECISION :: snu(jy),csu(jy),snv(jy),csv(jy)
  DOUBLE PRECISION :: cdfu(jy,-1:0),cdfv(jy,-1:0),ccou(jy),ccov(jy),cpgu(jy),cpgv
  DOUBLE PRECISION :: cst,ccn1(jy),ccn2(jy)
  DOUBLE PRECISION :: chvu1(jy),chv2(jy),chvu3(jy),chvu4(jy),chvu5(jy)
  DOUBLE PRECISION :: chvv1(jy),         chvv3(jy),chvv4(jy),chvv5(jy)
  DOUBLE PRECISION :: taux(ix,jy),tauy(ix,jy)
  INTEGER :: igrdh(ix,jy),igrdu(ix,jy),igrdv(ix,jy)
  INTEGER :: iu(ix*2,2),ju(ix*2,2),iv(jy*2,2),jv(jy*2,2),iumx(2),ivmx(2)
END MODULE commons

! **********************************************************************
FUNCTION fi2lon(i)
  USE commons
  IMPLICIT NONE
  INTEGER :: i
  DOUBLE PRECISION :: fi2lon
  fi2lon =  DBLE(i-2)*dx + x0
  RETURN
END FUNCTION fi2lon
! **********************************************************************
FUNCTION fj2lat(j)
  USE commons
  IMPLICIT NONE
  INTEGER :: j
  DOUBLE PRECISION :: fj2lat
  fj2lat = (DBLE(j-2)*dy+y0)*180.0d0/pi
  RETURN
END FUNCTION fj2lat
! **********************************************************************
FUNCTION lon2i(flon)
  USE commons
  IMPLICIT NONE
  INTEGER :: lon2i
  DOUBLE PRECISION :: flon
  lon2i = INT((flon-x0)/dx)+2
  RETURN
END FUNCTION lon2i
! **********************************************************************
FUNCTION lat2j(flat)
  USE commons
  IMPLICIT NONE
  INTEGER :: lat2j
  DOUBLE PRECISION :: flat
  lat2j = INT((flat-y0)/dy)+2
  RETURN
END FUNCTION lat2j
! **********************************************************************

! **********************************************************************
SUBROUTINE parameters
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  DOUBLE PRECISION :: g, torad

  ahmn  = 1.0d4         ! [m**2/s]
  ahmx  = 1.0d5         ! [m**2/s]
  bc    = 1.0d0         ! boundary condition (slippery=1.,viscous=-1.)
  dt1   = 900.0d0       ! [s]
  dx    = 1.0d0         ! [deg]
  dy    = 1.0d0         ! [deg]
  g     = 9.8d0         ! [m/s**2]
  hini  = 200.0d0       ! [m]
  re    = 6378000.0d0   ! [m]
  rho1  = 1026.0d0      ! [kg/m**3]
  rho2  = 1028.0d0      ! [kg/m**3]
  x0    = 100.0d0       ! [degree east]
  y0    = -50.0d0       ! [degree north]

  pi    = 4.0d0*ATAN(1.0d0)
  torad = pi/180.0d0
  dt2   = 2.0d0*dt1
  dx    = dx*torad
  dy    = dy*torad
  rg    = g*(rho2-rho1)/rho2
  omega = 2.0d0*pi/(24.0d0*3600.0d0)
  x0    = x0*torad
  y0    = y0*torad

  RETURN
END SUBROUTINE parameters

! **********************************************************************
SUBROUTINE initialization
  ! **********************************************************************
  USE commons
  IMPLICIT NONE

  u(:,:,:)   = 0.0d0
  v(:,:,:)   = 0.0d0
  h(:,:,:)   = 0.0d0
  taux(:,:)  = 0.0d0
  tauy(:,:)  = 0.0d0
  igrdh(:,:) = 0
  igrdu(:,:) = 0
  igrdv(:,:) = 0

  cdfu(:,:)= 0.0d0
  cdfv(:,:)= 0.0d0
  ccou(:)  = 0.0d0
  ccov(:)  = 0.0d0
  cpgu(:)  = 0.0d0
  ccn1(:)  = 0.0d0
  ccn2(:)  = 0.0d0
  chvu1(:) = 0.0d0
  chv2(:)  = 0.0d0
  chvu3(:) = 0.0d0
  chvu4(:) = 0.0d0
  chvu5(:) = 0.0d0
  chvv1(:) = 0.0d0
  chvv3(:) = 0.0d0
  chvv4(:) = 0.0d0
  chvv5(:) = 0.0d0

  iu(:,:) = 0
  ju(:,:) = 0
  iv(:,:) = 0
  jv(:,:) = 0

  RETURN
END SUBROUTINE initialization

! **********************************************************************
SUBROUTINE coefficients
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: j

  DOUBLE PRECISION :: dmdfu(jy), dmdfv(jy), fj2lat
  ! ------------------------------ trigonometric functions
  DO j=1,jy
     snu(j) = SIN(y0+DBLE(j-2)*dy)
     snv(j) = SIN(y0+DBLE(j-2)*dy-0.5d0*dy)
     csu(j) = COS(y0+DBLE(j-2)*dy)
     csv(j) = COS(y0+DBLE(j-2)*dy-0.5d0*dy)
  ENDDO
  ! ------------------------------ ah
  DO j=1,jy
     IF ( ABS(fj2lat(j))>=40.0d0 ) THEN
        ah(j) = ahmn+(ahmx-ahmn)*( ABS(fj2lat(j))-40.0d0)/10.0d0
     ELSE
        ah(j) = ahmn
     ENDIF
     !     PRINT *, "j=", j, " lat=", fj2lat(j), " ah=", ah(j)
  ENDDO
  ! ------------------------------ coefficients for dufort-frankel scheme
  ! ------------------------------                                 part1
  DO j=1,jy
     dmdfu(j) = ah(j) / ((re*csu(j)*dx)**2.0d0)+ah(j)/((re*dy)**2.0d0)
     dmdfv(j) = ah(j) / ((re*csv(j)*dx)**2.0d0)+ah(j)/((re*dy)**2.0d0)
  ENDDO
  ! ------------------------------ coefficients for dufort-frankel scheme
  ! ------------------------------                                 part2
  DO j=1,jy
     cdfu(j,-1) = 1.0d0 / (1.0d0+dt2*dmdfu(j))
     cdfu(j, 0) = 1.0d0 / (1.0d0+dt1*dmdfu(j))
     cdfv(j,-1) = 1.0d0 / (1.0d0+dt2*dmdfv(j))
     cdfv(j, 0) = 1.0d0 / (1.0d0+dt1*dmdfv(j))
  ENDDO
  ! ------------------------------ coefficients for colioris term
  DO j=1,jy
     ccou(j) =  omega*snu(j)*0.5d0
     ccov(j) = -omega*snv(j)*0.5d0
  ENDDO
  ! ------------------------------ coefficients for pressure gradient term
  DO j=1,jy
     cpgu(j) = -rg / (re*csu(j)*dx)
  ENDDO
  cpgv = -rg / (re*dy)
  ! ------------------------------ coefficient for stress term
  cst = 0.5d0 / rho1
  ! ------------------------------ coefficients for
  ! ------------------------------     horizontal eddy viscosity term
  DO j=1,jy
     chvu1(j) =  ah(j) / ((re*csu(j)*dx)**2.0d0)
     chv2(j)  =  ah(j) / ((re*dy)**2.0d0)
     chvu3(j) = -(ah(j)*snu(j)) / (re**2.0d0*csu(j)*2.0d0*dy)
     chvu4(j) = -(ah(j)*snu(j)) / (re**2.0d0*csu(j)**2.0d0*dx)
     chvu5(j) = -ah(j) / ((re*csu(j))**2.0d0)
     chvv1(j) =  ah(j) / ((re*csv(j)*dx)**2)
     chvv3(j) = -(ah(j)*snv(j)) / (re**2.0d0*csv(j)*2.0d0*dy)
     chvv4(j) =  (ah(j)*snv(j)) / (re**2.0d0*csv(j)**2.0d0*dx)
     chvv5(j) = -ah(j) / ((re*csv(j))**2.0d0)
  ENDDO
  ! ------------------------------ coefficients for continuity equation
  DO j=1,jy
     ccn1(j) = 1.0d0 / (re*csu(j)*dx)
     ccn2(j) = 1.0d0 / (re*csu(j)*dy)
  ENDDO
  !
  RETURN
END SUBROUTINE coefficients

! **********************************************************************
SUBROUTINE mapping
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j

  INTEGER :: land(ix,jy)
  ! ------------------------------ reading bottom topography
  ! ------------------------------ land-sea map(0=sea,1=land)
  DO j=1,jy
     DO i=1,ix
        land(i,j) = 1
     ENDDO
  ENDDO
  !
  DO j=1,jy
     !!  do j=2,jy-1
     DO i=2,ix-1
        land(i,j) = 0
     ENDDO
  ENDDO
  !
  !     open(50,file='../lib/poem.map',status='old')
  !     do j=jy-1,2,-1
  !        read(50,'(171i1)') (land(i,j),i=2,ix-1)
  !     enddo
  !     close(50)
  ! ------------------------------ setting calculation grid point for h,u,v(1=calc)
  ! ------------------------------ and setting initial upper layer thickness
  DO j=1,jy
     DO i=2,ix-1
        IF ( land(i,j)==0 ) THEN
           h(i,j, 0) = hini
           h(i,j,-1) = hini
           h(i,j,-2) = hini
           igrdh(i,j) = 1
        ENDIF
        IF ( land(i-1,j)==0 .AND. land(i,j)==0 ) igrdu(i,j) = 1
        IF ( land(i,j-1)==0 .AND. land(i,j)==0 ) igrdv(i,j) = 1
     ENDDO
  ENDDO
  ! ------------------------------ checking boundary part1
  DO j=2,jy
     DO i=2,ix
        IF (igrdh(i-1,j)==0 .AND. igrdh(i,j)==1) igrdu(i,j) = 0
        IF (igrdh(i-1,j)==1 .AND. igrdh(i,j)==0) igrdu(i,j) = 0
        IF (igrdh(i,j-1)==0 .AND. igrdh(i,j)==1) igrdv(i,j) = 0
        IF (igrdh(i,j-1)==1 .AND. igrdh(i,j)==0) igrdv(i,j) = 0
     ENDDO
  ENDDO
  ! ------------------------------ checking boundary part2
  iumx(1) = 0
  iumx(2) = 0
  ivmx(1) = 0
  ivmx(2) = 0

  DO j=2,jy-1
     DO i=2,ix-1
        IF (igrdh(i,j)==1) THEN

           ! -------------------- about u with south wall
           IF (igrdh(i-1,j-1)==0 .AND. igrdh(i,j-1)==0) THEN
              iumx(1)=iumx(1)+1
              iu(iumx(1),1)=i
              ju(iumx(1),1)=j
           ENDIF
           ! -------------------- about u with north wall
           IF (igrdh(i-1,j+1)==0 .AND. igrdh(i,j+1)==0) THEN
              iumx(2)=iumx(2)+1
              iu(iumx(2),2)=i
              ju(iumx(2),2)=j
           ENDIF
           ! -------------------- about v with west wall
           IF (igrdh(i-1,j-1)==0 .AND. igrdh(i-1,j)==0) THEN
              ivmx(1)=ivmx(1)+1
              iv(ivmx(1),1)=i
              jv(ivmx(1),1)=j
           ENDIF
           ! -------------------- about v with east wall
           IF (igrdh(i+1,j-1)==0 .AND. igrdh(i+1,j)==0) THEN
              ivmx(2)=ivmx(2)+1
              iv(ivmx(2),2)=i
              jv(ivmx(2),2)=j
           ENDIF

        ENDIF
     ENDDO
  ENDDO
  ! ------------------------------ Validation
  DO i=1,ix
     IF (igrdh(i, 1)>0 ) &
          &     PRINT *, 'Ocean(h) on the South Wall!', 'i=', i, 'j=', j
     IF (igrdh(i,jy)>0 ) &
          &     PRINT *, 'Ocean(h) on the North Wall!', 'i=', i, 'j=', j
  ENDDO

  DO j = 1, jy
     IF (igrdh( 1,j)>0 ) &
          &     PRINT *, 'Ocean(h) on the West Wall!', 'i=', i, 'j=', j
     IF (igrdh(ix,j)>0 ) &
          &     PRINT *, 'Ocean(h) on the East Wall!', 'i=', i, 'j=', j
  ENDDO
  !
  DO j=2,jy-1
     DO i=2,ix-1
        !
        IF (igrdh(i-1,j)==0 .AND. igrdh(i,j)==0 .AND. igrdh(i,j)==1) &
             &       PRINT *, 'Ocean(u) in the Land!', 'i=', i, 'j=', j
        IF (igrdh(i-1,j)==1 .AND. igrdh(i,j)==1 .AND. igrdu(i,j)==0) &
             &       PRINT *, 'Land(u) in the Ocean!', 'i=', i, 'j=', j
        !
        IF (igrdh(i,j-1)==0 .AND. igrdh(i,j)==0 .AND. igrdv(i,j)==1) &
             &       PRINT *, 'Ocean(v) in the Land!', 'i=', i, 'j=', j
        IF (igrdh(i,j-1)==1 .AND. igrdh(i,j)==1 .AND. igrdv(i,j)==0) &
             &       PRINT *, 'Land(v) in the Ocean!', 'i=', i, 'j=', j
        !
        IF (igrdu(i,j)==0 .AND. igrdu(i+1,j)==0 .AND. igrdh(i,j)==1) &
             &       PRINT *, 'Ocean(h) in the Land!', 'i=', i, 'j=', j
        IF (igrdu(i,j)==1 .AND. igrdu(i+1,j)==1 .AND. igrdh(i,j)==0) &
             &       PRINT *, 'Land(h) in the Ocean!', 'i=', i, 'j=', j
        !
        IF (igrdv(i,j)==0 .AND. igrdv(i,j+1)==0 .AND. igrdh(i,j)==1) &
             &       PRINT *, 'Ocean(h) in the Land!', 'i=', i, 'j=', j
        IF (igrdv(i,j)==1 .AND. igrdv(i,j+1)==1 .AND. igrdh(i,j)==0)  &
             &       PRINT *, 'Land(h) in the Ocean!', 'i=', i, 'j=', j
        !
     ENDDO
  ENDDO
  !
  OPEN(61,file='../dat/land.map')
  DO j=jy,1,-1
     WRITE (61,'(203i1)') (land(i,j),i=1,ix)
  ENDDO
  CLOSE(61)
  OPEN(62,file='../dat/igrdu.map')
  DO j=jy,1,-1
     WRITE (62,'(203i1)') (igrdu(i,j),i=1,ix)
  ENDDO
  CLOSE(62)
  OPEN(63,file='../dat/igrdv.map')
  DO j=jy,1,-1
     WRITE (63,'(203i1)') (igrdv(i,j),i=1,ix)
  ENDDO
  CLOSE(63)
  OPEN(64,file='../dat/igrdh.map')
  DO j=jy,1,-1
     WRITE (64,'(203i1)') (igrdh(i,j),i=1,ix)
  ENDDO
  CLOSE(64)

  !
  RETURN
END SUBROUTINE mapping

! **********************************************************************
SUBROUTINE forcing
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j
  DOUBLE PRECISION :: r

  ! ------------------------------ wind stress
  ! ------------------------------ reading bottom topography
  !      open(51,file='coads-wndstr.dat',
  !     +        form='unformatted',
  !     +        access='direct',
  !     +        recl=(ix-2)*(jy-2)*4)
  !      read(51,rec=1) ((taux(i,j),i=2,ix-1),j=2,jy-1)
  !      read(51,rec=2) ((tauy(i,j),i=2,ix-1),j=2,jy-1)
  !      close(51)
  !
  DO j=1,jy
     DO i=1,ix
        !     taux(i,j)=  0.01d0*cos(2.0d0*pi*dble(j-11)/dble(103-1)) ! 2gyres
        !      taux(i,j)=  0.01d0*COS(2.0d0*pi*DBLE(j-1)/DBLE(jy-1)) ! 2gyres
        !     taux(i,j)= -0.01d0*cos(pi*dble(j-1)/dble(jy-1))      ! 1gyre
        !     taux(i,j)= -0.01d0
        !taux(i,j)=  0.00d0

        tauy(i,j)=  0.00d0
     ENDDO
  ENDDO
  ! ------------------------------ big water column at center of the ocean
  DO j=INT(jy/2)-10,INT(jy/2)+10
     DO i=INT(ix/2)-10,INT(ix/2)+10
        r = SQRT((DBLE(ix/2+1-i))**2+(DBLE(jy/2+1-j))**2)
        IF ( r<=10.0d0 ) THEN
           h(i,j,0) = h(i,j,0)+20.0d0
           h(i,j,-1) = h(i,j,-1)+20.0d0
           h(i,j,-2) = h(i,j,-2)+20.0d0
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE forcing

! **********************************************************************
SUBROUTINE calculation(it)
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j, n, it
  DOUBLE PRECISION :: a, dt, cl

  !
  IF (MOD(it,100)==1) THEN
     dt = dt1
     n  = 0
  ELSE
     dt = dt2
     n  = -1
  ENDIF
  ! ------------------------------ giving time-variation to wind stress
  a = 1.0d0
  IF (it<=INT(10.0d0*24.0d0*3600.0d0/dt1)) &
       &   a = DBLE(it)*dt1/(10.0d0*24.0d0*3600.0d0)
  ! ------------------------------ calculating momentum equation
  DO j=2,jy-1
     DO i=2,ix-1
        !
        u(i,j,1) = cdfu(j,n)*(u(i,j,n)+dt*( &
             &     ccou(j)*(v(i-1,j,0)+v(i,j,0)+v(i-1,j+1,0)+v(i,j+1,0)) &
             &    +cpgu(j)*0.5d0*(h(i,j,0)+h(i-1,j,0))*(h(i,j,0)-h(i-1,j,0)) &
             &    +a*cst*(taux(i-1,j)+taux(i,j)) &
             &    +chvu1(j)*(u(i+1,j,0)+u(i-1,j,0)-u(i,j,-1)) &
             &    +chv2(j)*(u(i,j+1,0)+u(i,j-1,0)-u(i,j,-1)) &
             &    +chvu3(j)*(u(i,j+1,n)-u(i,j-1,n)) &
             &    +chvu4(j)*(v(i,j+1,n)+v(i,j,n)-v(i-1,j+1,n)-v(i-1,j,n)) &
             &    +chvu5(j)*u(i,j,n) &
             & ))
        !
        v(i,j,1) = cdfv(j,n)*(v(i,j,n)+dt*( &
             &     ccov(j)*(u(i,j-1,0)+u(i+1,j-1,0)+u(i,j,0)+u(i+1,j,0)) &
             &    +cpgv*0.5d0*(h(i,j,0)+h(i,j-1,0))*(h(i,j,0)-h(i,j-1,0)) &
             &    +a*cst*(tauy(i,j-1)+tauy(i,j)) &
             &    +chvv1(j)*(v(i+1,j,0)+v(i-1,j,0)-v(i,j,-1)) &
             &    +chv2(j)*(v(i,j+1,0)+v(i,j-1,0)-v(i,j,-1)) &
             &    +chvv3(j)*(v(i,j+1,n)-v(i,j-1,n)) &
             &    +chvv4(j)*(u(i+1,j,n)+u(i+1,j-1,n)-u(i,j,n)-u(i,j-1,n)) &
             &    +chvv5(j)*v(i,j,n) &
             & ))
        !
        h(i,j,1) = h(i,j,n)-dt*( &
             &     ccn1(j)*(u(i+1,j,0)-u(i,j,0)) &
             &    +ccn2(j)*(v(i,j+1,0)*csv(j+1)-v(i,j,0)*csv(j)) &
             & )
     ENDDO
  ENDDO
  !
  ! ------------------------------ Apply Orlanski radiation condition
  ! ------------------------------ at Northern open boundary
  DO i=1,ix
     IF (igrdu(i,jy)==0) CYCLE
     cl = ( u(i,jy-1,-2)-u(i,jy-1,0) ) / ( u(i,jy-1,0)+u(i,jy-1,-2)-2.0d0*u(i,jy-2,-1) )
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        u(i,jy,1) = u(i,jy,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        u(i,jy,1) = u(i,jy-1,0)
     ELSE                ! 0 < μ < 1
        u(i,jy,1) = ( u(i,jy,-1)*(1.0d0-cl)+2.0d0*cl*u(i,jy-1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO
  !
  DO i=1,ix
     IF (igrdv(i,jy)==0) CYCLE
     cl = ( v(i,jy-1,-2)-v(i,jy-1,0) ) / ( v(i,jy-1,0)+v(i,jy-1,-2)-2.0d0*v(i,jy-2,-1) )
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        v(i,jy,1) = v(i,jy,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        v(i,jy,1) = v(i,jy-1,0)
     ELSE                ! 0 < μ < 1
        v(i,jy,1) = ( v(i,jy,-1)*(1.0d0-cl)+2.0d0*cl*v(i,jy-1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO
  !
  DO i=1,ix
     IF (igrdh(i,jy)==0) CYCLE
     cl = ( h(i,jy-1,-2)-h(i,jy-1,0) ) / ( h(i,jy-1,0)+h(i,jy-1,-2)-2.0d0*h(i,jy-2,-1) )
     !
     !     IF (i==102) PRINT *, i, cl, h(i,jy-1,-2), h(i,jy-1,0), h(i,jy-1,0), h(i,jy-1,-2), h(i,jy-2,-1)
     !
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        h(i,jy,1) = h(i,jy,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        h(i,jy,1) = h(i,jy-1,0)
     ELSE                ! 0 < μ < 1
        h(i,jy,1) = ( h(i,jy,-1)*(1.0d0-cl)+2.0d0*cl*h(i,jy-1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO

  ! ------------------------------ Apply Orlanski radiation condition
  ! ------------------------------ at Southern open boundary
  DO i=1,ix
     IF (igrdu(i,1)==0) CYCLE
     cl = ( u(i,1+1,-2)-u(i,1+1,0) ) / ( u(i,1+1,0)+u(i,1+1,-2)-2.0d0*u(i,1+2,-1) )
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        u(i,1,1) = u(i,1,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        u(i,1,1) = u(i,1+1,0)
     ELSE                ! 0 < μ < 1
        u(i,1,1) = ( u(i,1,-1)*(1.0d0-cl)+2.0d0*cl*u(i,1+1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO
  !
  DO i=1,ix
     IF (igrdv(i,1)==0) CYCLE
     cl = ( v(i,1+1,-2)-v(i,1+1,0) ) / ( v(i,1+1,0)+v(i,1+1,-2)-2.0d0*v(i,1+2,-1) )
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        v(i,1,1) = v(i,1,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        v(i,1,1) = v(i,1+1,0)
     ELSE                ! 0 < μ < 1
        v(i,1,1) = ( v(i,1,-1)*(1.0d0-cl)+2.0d0*cl*v(i,1+1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO
  !
  DO i=1,ix
     IF (igrdh(i,1)==0) CYCLE
     cl = ( h(i,1+1,-2)-h(i,1+1,0) ) / ( h(i,1+1,0)+h(i,1+1,-2)-2.0d0*h(i,1+2,-1) )
     IF (cl<0 .OR. isnan(cl) ) THEN      ! μ = 0
        h(i,1,1) = h(i,1,-1)
     ELSE IF (cl>1) THEN ! μ = 1
        h(i,1,1) = h(i,1+1,0)
     ELSE                ! 0 < μ < 1
        h(i,1,1) = ( h(i,1,-1)*(1.0d0-cl)+2.0d0*cl*h(i,1+1,0) ) / ( 1.0d0+cl )
     ENDIF
  ENDDO
  !
  RETURN
END SUBROUTINE calculation

! **********************************************************************
SUBROUTINE energy(it)
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j, n, it
  DOUBLE PRECISION :: eng, have

  !
  eng=0.0d0
  have=0.0d0
  !
  DO n=1,iumx(1)
     eng = eng+0.5d0*(u(iu(n,1),ju(n,1),0)/ &
          &              (0.5d0*(h(iu(n,1)-1,ju(n,1),0)+h(iu(n,1),ju(n,1),0))))**2.0d0
  ENDDO
  DO n=1,iumx(2)
     eng = eng+0.5d0*(u(iu(n,2),ju(n,2),0)/ &
          &              (0.5d0*(h(iu(n,2)-1,ju(n,2),0)+h(iu(n,2),ju(n,2),0))))**2.0d0
  ENDDO
  DO n=1,ivmx(1)
     eng = eng+0.5d0*(v(iv(n,1),jv(n,1),0)/ &
          &              (0.5d0*(h(iv(n,1),jv(n,1)-1,0)+h(iv(n,1),jv(n,1),0))))**2.0d0
  ENDDO
  DO n=1,ivmx(2)
     eng = eng+0.5d0*(v(iv(n,2),jv(n,2),0)/ &
          &              (0.5d0*(h(iv(n,2),jv(n,2)-1,0)+h(iv(n,2),jv(n,2),0))))**2.0d0
  ENDDO
  !
  DO j=2,jy-1
     DO i=2,ix-1
        IF (igrdh(i,j)/=0) eng = eng+DBLE(igrdh(i,j)) &
             &                *0.5d0 *((0.5d0*(u(i,j,0)+u(i+1,j,0))/h(i,j,0))**2.0d0 &
             &                        +(0.5d0*(v(i,j,0)+v(i,j+1,0))/h(i,j,0))**2.0d0)
     ENDDO
  ENDDO
  !
  PRINT *, DBLE(it)/(24.0d0*3600.0d0/dt1),eng,SUM(h(:,:,0))/SUM(igrdh)
  !!  print '(7e10.3)', dble(it)/(24.0d0*3600.0d0/dt1),v(182,101,0),h(182,101,0),v(182,102,0),h(182,102,0),v(182,103,0),h(182,103,0)
  !
  RETURN
END SUBROUTINE energy

! **********************************************************************
SUBROUTINE nextstep()
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j, n

  !
  DO j=1,jy
     DO i=1,ix
        u(i,j,-2) = u(i,j,-1)*DBLE(igrdu(i,j))
        u(i,j,-1) = u(i,j, 0)*DBLE(igrdu(i,j))
        u(i,j, 0) = u(i,j, 1)*DBLE(igrdu(i,j))

        v(i,j,-2) = v(i,j,-1)*DBLE(igrdv(i,j))
        v(i,j,-1) = v(i,j, 0)*DBLE(igrdv(i,j))
        v(i,j, 0) = v(i,j, 1)*DBLE(igrdv(i,j))

        h(i,j,-2) = h(i,j,-1)*DBLE(igrdh(i,j))
        h(i,j,-1) = h(i,j, 0)*DBLE(igrdh(i,j))
        h(i,j, 0) = h(i,j, 1)*DBLE(igrdh(i,j))
     ENDDO
  ENDDO
  ! ------------------------------ southern boundary
  DO n=1,iumx(1)
     u(iu(n,1),ju(n,1)-1, 0) = u(iu(n,1),ju(n,1), 0)*bc
     u(iu(n,1),ju(n,1)-1,-1) = u(iu(n,1),ju(n,1),-1)*bc
     u(iu(n,1),ju(n,1)-1,-2) = u(iu(n,1),ju(n,1),-2)*bc
  ENDDO
  ! ------------------------------ northern boundary
  DO n=1,iumx(2)
     u(iu(n,2),ju(n,2)+1, 0) = u(iu(n,2),ju(n,2), 0)*bc
     u(iu(n,2),ju(n,2)+1,-1) = u(iu(n,2),ju(n,2),-1)*bc
     u(iu(n,2),ju(n,2)+1,-2) = u(iu(n,2),ju(n,2),-2)*bc
  ENDDO
  ! ------------------------------ western boundary
  DO n=1,ivmx(1)
     v(iv(n,1)-1,jv(n,1), 0) = v(iv(n,1),jv(n,1), 0)*bc
     v(iv(n,1)-1,jv(n,1),-1) = v(iv(n,1),jv(n,1),-1)*bc
     v(iv(n,1)-1,jv(n,1),-2) = v(iv(n,1),jv(n,1),-2)*bc
  ENDDO
  ! ------------------------------ eastern boundary
  DO n=1,ivmx(2)
     v(iv(n,2)+1,jv(n,2), 0) = v(iv(n,2),jv(n,2), 0)*bc
     v(iv(n,2)+1,jv(n,2),-1) = v(iv(n,2),jv(n,2),-1)*bc
     v(iv(n,2)+1,jv(n,2),-2) = v(iv(n,2),jv(n,2),-2)*bc
  ENDDO
  ! ------------------------------ open boundary at northernmost & southernmost boundary
  !  DO i=1,ix
  !
  !     u(i, 1, 0) = u(i, 2  , 0)
  !     u(i, 1,-1) = u(i, 2  ,-1)
  !     u(i, 1,-2) = u(i, 2  ,-2)
  !
  !     u(i,jy, 0) = u(i,jy-1, 0)
  !     u(i,jy,-1) = u(i,jy-1,-1)
  !     u(i,jy,-2) = u(i,jy-1,-2)
  !
  !     v(i, 1, 0) = v(i, 3  , 0)
  !     v(i, 1,-1) = v(i, 3  ,-1)
  !     v(i, 1,-2) = v(i, 3  ,-2)
  !
  !     v(i, 2, 0) = v(i, 3  , 0)
  !     v(i, 2,-1) = v(i, 3  ,-1)
  !     v(i, 2,-2) = v(i, 3  ,-2)
  !
  !     v(i,jy, 0) = v(i,jy-1, 0)
  !     v(i,jy,-1) = v(i,jy-1,-1)
  !     v(i,jy,-2) = v(i,jy-1,-2)
  !
  !     h(i, 1, 0) = h(i, 2  , 0)
  !     h(i, 1,-1) = h(i, 2  ,-1)
  !     h(i, 1,-2) = h(i, 2  ,-2)
  !
  !     h(i,jy, 0) = h(i,jy-1, 0)
  !     h(i,jy,-1) = h(i,jy-1,-1)
  !     h(i,jy,-2) = h(i,jy-1,-2)
  !  ENDDO
  !
  RETURN
END SUBROUTINE nextstep

! **********************************************************************
SUBROUTINE information(iend,iout)
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: iend, iout

  !
  PRINT '("*********************")'
  PRINT '("* model information *")'
  PRINT '("*********************")'
  PRINT '("* -----end of calculations & output data interval")'
  PRINT '("the model was started from a state of rest.")'
  PRINT '("end of calculations  =",i10,"steps ",f8.2,"days")', iend,DBLE(iend)*dt1/(24.0d0*3600.0d0)
  PRINT '("output interval      =",i10,"steps ",f8.2,"days")', iout,DBLE(iout)*dt1/(24.0d0*3600.0d0)
  PRINT '("* -----parameters")'
  PRINT '("dt                   =",f6.1,"(s)")', dt1
  PRINT '("2*dt                 =",f6.1,"(s)")', dt2
  PRINT '("dx                   =",e10.3,"(rad)")', dx
  PRINT '("dy                   =",e10.3,"(rad)")', dy
  PRINT '("ah(min)              =",e10.3,"(m**2/sec)")', ahmn
  PRINT '("ah(max)              =",e10.3,"(m**2/sec)")', ahmx
  PRINT '("reduced gravity      =",e10.3,"(m/s**2)")', rg
  PRINT '("omega                =",e10.3,"(1/s)")', omega
  PRINT '("radius of the earth  =",e10.3,"(m)")', re
  PRINT '("rho1                 =",f6.1,"(kg/m**3)")', rho1
  PRINT '("rho2                 =",f6.1,"(kg/m**3)")', rho2
  PRINT '("initial layer thickness =",f6.1,"(m)")', hini
  PRINT '("* -----boundary condition")'
  PRINT '("(1.0=slippery,-1.0=viscous)...",f4.1)', bc
  PRINT '("* -----number of boundary points")'
  PRINT '("* -----   which needs data replacing")'
  PRINT '("on u-point boundary at south =",i4)', iumx(1)
  PRINT '("           boundary at north =",i4)', iumx(2)
  PRINT '("on v-point boundary at west  =",i3)', ivmx(1)
  PRINT '("           boundary at east  =",i3)', ivmx(2)
  PRINT '("* --------------------end of informations")'
  !
  RETURN
END SUBROUTINE information

! **********************************************************************
SUBROUTINE output(irec)
  ! **********************************************************************
  USE commons
  IMPLICIT NONE
  INTEGER :: i, j, irec
  REAL :: uu(ix,jy), vv(ix,jy)
  !
  OPEN(60,file='../dat/output.dat', &
       &         form='unformatted', access='direct', recl=(ix-2)*(jy-2)*4)
  !
  DO j=2,jy-1
     DO i=2,ix-1
        IF (igrdh(i,j)/=0) THEN
           uu(i,j) = REAL(igrdh(i,j))*0.5*REAL((u(i,j,0)+u(i+1,j,0))/h(i,j,0))
           vv(i,j) = REAL(igrdh(i,j))*0.5*REAL((v(i,j,0)+v(i,j+1,0))/h(i,j,0))
        ELSE
           uu(i,j) = 0.0
           vv(i,j) = 0.0
        ENDIF
     ENDDO
  ENDDO
  !
  ! print *, 'irec->',irec
  WRITE(60,rec=irec)   ((uu(i,j),i=2,ix-1),j=2,jy-1)
  WRITE(60,rec=irec+1) ((vv(i,j),i=2,ix-1),j=2,jy-1)
  WRITE(60,rec=irec+2) ((REAL(h(i,j,0)),i=2,ix-1),j=2,jy-1)
  irec=irec+3
  !
  CLOSE(60)
  !
  RETURN
END SUBROUTINE output

! **********************************************************************
! **********************************************************************
PROGRAM yum
  ! **********************************************************************
  ! **********************************************************************

  USE commons
  IMPLICIT NONE
  INTEGER :: it, iend, iout, irec

  ! ---------------------------- define model run
  iend = 10000          ! end of run [days]
  iout = 10            ! output interval [days]
  ! ---------------------------- initialize
  CALL initialization
  ! ---------------------------- parameters
  CALL parameters

  iend = iend *INT(24.0d0*3600.0d0/dt1)
  iout = iout *INT(24.0d0*3600.0d0/dt1)
  irec = 1
  ! ---------------------------- coefficients
  CALL coefficients
  ! ---------------------------- mapping
  CALL mapping
  ! ---------------------------- bottom & wind
  CALL forcing
  ! ---------------------------- output model information
  CALL information(iend,iout)

  ! ------------------ output initial state & calculating momentum energy
  CALL energy(it)
  CALL output(irec)

  ! ---------------------------- loop start

  DO it=1,iend

     ! ------------------ calculating momentum equations & continuity equations
     CALL calculation(it)
     ! ------------------ output result & calculating momentum energy
     IF (MOD(it,iout)==0) THEN
        CALL energy(it)
        CALL output(irec)
     ENDIF
     ! ------------------ for next step
     CALL nextstep()

     ! ---------------------------- loop end
  ENDDO

  ! —————————————— end
END PROGRAM yum


! **********************************************************************
! Grid Position
!
!     h   u - h   u - h    j+1
!     |     \ |     \ |
!     v       v       v
!
!     h   u - h   u - h    j
!     |     \ |     \ |
!     v       v       v
!
!     h   u - h   u - h    j-1
!
!    i-1      i      i+1
!
! **********************************************************************
! Variables & Parameters
!     u(i,j,k)......transport of eastward component (k=1:past,k=2:now,k=3:future)
!     v(i,j,k)......transport of eastward component (k=1:past,k=2:now,k=3:future)
!     h(i,j,k)......upper layer thickness(k=1:past,k=2:now,k=3:future)
!     dt............time interval (s)
!     dx............grid interval of eastward component (deg)
!     dy............grid interval of northward component (deg)
!     omega.........angular rotation of the earth (1/s)
!     re............radius of the earth (m)
!     ah(j).........horizontal eddy viscosity coefficient (m**2/s**2)
!     fu(j),fv(j)...colioris parameter
!     h(i,j)........upper layer thickness (m)
!     rho1..........density of upper layer (kg/m**3)
!     rho2..........density of lower layer (kg/m**3)
!     g.............acceleration of gravity=(9.8m/s**2)
!     rg............reduced gravity  (m/s**2)
!     taux(i,j).....wind stress of eastward component (m**2/s)
!     tauy(i,j).....wind stress of northward component (m**2/s)
!     cl............Courant Number　for Orlanski Radiation Condition
!
! **********************************************************************
! Orlanski Radiation Condition
!                    n-2    n
!                   φ    - φ
!                    B∓1    B∓1
!     μ = CL = ----------------------
!                n      n-2      n-1
!               φ    + φ    - 2•φ
!                B∓1    B∓1      B∓2
!
!     μ=CL if 0<CL<1:
!                 n-1                n
!                φ   •( 1-μ ) + 2•μ•φ
!         n+1     B                  B∓1
!        φ    = -------------------------
!         B               1 + μ
!
!     μ=0 if CL<0:
!         n+1     n-1
!        φ    =  φ
!         B       B
!
!     μ=1 if  CL>1:
!         n+1     n
!        φ    =  φ
!         B       B∓1
!
! **********************************************************************
