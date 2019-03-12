!===============================================================================
!===============================================================================
MODULE radau_mod
  implicit none
  private
  type, public:: radau_t
  contains
    procedure, nopass:: init
  end type
  type(radau_t), public:: radau
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (k,xi,wi)
!
!  Routine to return weights and abscissa points for Gauss-Radau
!  integration in the interval [a; b]. The transformation from
!  the inherent interval of [-1; 1] is performed via c and d.
!  The Radau integration differs from a normal Gauss integration,
!  in that one of the end-points are included.
!  The number of integration points, k, can be chosen between
!  1 and 10, and the order of the interpolating function is 2*k-1.
!
!  Update history:
!  ---------------
!  29.08.2006 Coded/RT
!  31.08.2006 Fixed single/double precision ambiguity and changed from
!             x=0 to x=1 being included as a node/RT
!  15.01.2007 Ported from F77 to F90/RT
!
!  Known bugs: Hardcoded whether to use upper (coded) or lower end-point
!              in the integration.
!
  implicit none
  integer, parameter:: kmax = 10
  integer::              k
  real, dimension(k)::  wi, xi
  real::                 a, b
  real(kind=8), dimension(55)::  xn = (/ &
          -1.00000000000000d0,                                      &!  1
          -1.00000000000000d0, .33333333333333d0,                   &!  2
          -1.00000000000000d0,-.28989794855663d0, .68989794855663d0,&!  3
          -1.00000000000000d0,-.57531892352169d0, .18106627111853d0,&!  4
            .82282408097459d0,                                      &
          -1.00000000000000d0,-.72048027131244d0,-.16718086473783d0,&!  5
            .44631397272376d0, .88579160777096d0,                   &
          -1.00000000000000d0,-.80292982840235d0,-.39092854670727d0,&!  6
            .12405037950522d0, .60397316425279d0, .92038028589707d0,&
          -1.00000000000000d0,-.85389134263949d0,-.53846772406011d0,&!  7
           -.11734303754309d0, .32603061943769d0, .70384280066303d0,&
            .94136714568042d0,                                      &
          -1.00000000000000d0,-.88747487892616d0,-.63951861652622d0,&!  8
           -.29475056577366d0, .09430725266111d0, .46842035443083d0,&
            .77064189367819d0, .95504122712258d0,                   &
          -1.00000000000000d0,-.91073208942006d0,-.71126748591571d0,&!  9
           -.42635048571114d0,-.09037336960685d0, .25613567083346d0,&
            .57138304120873d0, .81735278420041d0, .96444016970528d0,&
          -1.00000000000000d0,-.92748437423359d0,-.76384204242000d0,&! 10
           -.52564603037008d0,-.23623446939059d0, .07605919783798d0,&
            .38066484014473d0, .64776668767401d0, .85122522058160d0,&
            .97117518070224d0/)
  real(kind=8), dimension(55)::  wn = (/ &
           2.00000000000000d0,                                      &!  1
           0.50000000000000d0,1.50000000000005d0,                   &!  2
           0.22222222222222d0,1.02497165237682d0,0.75280612540100d0,&!  3
           0.12500000000000d0,0.65768863996010d0,0.77638693768637d0,&!  4
           0.44092442235367d0,                                      &
           0.08000000000000d0,0.44620780216715d0,0.62365304595145d0,&!  5
           0.56271203029887d0,0.28742712158247d0,                   &
           0.05555555555556d0,0.31964075322051d0,0.48538718846895d0,&!  6
           0.52092678318961d0,0.41690133431186d0,0.20158838525329d0,&
           0.04081632653061d0,0.23922748922532d0,0.38094987364422d0,&!  7
           0.44710982901453d0,0.42470377900601d0,0.31820423146727d0,&
           0.14898847111224d0,                                      &
           0.03125000000000d0,0.18535815480298d0,0.30413062064680d0,&!  8
           0.37651754538912d0,0.39157216745253d0,0.34701479563444d0,&
           0.24964790132996d0,0.11450881474417d0,                   &
           0.02469135802469d0,0.14765401904631d0,0.24718937820459d0,&!  9
           0.31684377567046d0,0.34827300277294d0,0.33769396697589d0,&
           0.28638669635730d0,0.20055329802460d0,0.09071450492309d0,&
           0.02000000000000d0,0.12029667055749d0,0.20427013187899d0,&! 10
           0.26819483784117d0,0.30585928772442d0,0.31358245722693d0,&
           0.29061016483286d0,0.23919343171435d0,0.16437601273700d0,&
           0.07361700548689d0/)
  integer, dimension(kmax+1):: k0 = (/1,2,4,7,11,16,22,29,37,46,56/)
  real c,d
  integer i
  
  a = 0.0
  b = 1.0
  if (k.lt.1) then
    print *,' WARNING(radaui): k=',k,'< 1 does not make sense!'
    print *,'                  Setting k=1.'
    k = 1
  endif
  if (k.gt.kmax) then
    print *,' WARNING(radaui): k=',k,'>',kmax,' is not implemented.'
    print *,'                  Setting k=',kmax,'.'
    k = kmax
  endif
  c = .5*(b+a)
  d = .5*(b-a)
!
! Loop through the number of points in the integration
  do i=1, k
!
! Using this form instead of [c + d*xn(k0(k)+i-1)] makes
! x=1 (or b) the predefined node instead of x=-1 (or a).
    xi(i) = 1.0 - (c + d*xn(k0(k+1)-i))
    wi(i) =            d*wn(k0(k+1)-i)
  enddo

  return
END SUBROUTINE

END MODULE radau_mod
