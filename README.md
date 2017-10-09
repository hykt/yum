# yum
Yokota's Undeveloped Model

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

