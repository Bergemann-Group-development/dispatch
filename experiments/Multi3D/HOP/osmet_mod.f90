MODULE osmet_mod
    USE os_mod
    USE prec_mod
    USE op_math_mod
    implicit none

    real(dp), parameter :: ln10 = log(10.0)
    integer, parameter :: mnu=154000, mtem=17, mp6=18

    type :: osmet_type
        real(dp) :: nu(mnu)
        real(dp), dimension(:), allocatable :: p6, temp
        real(4), dimension(:,:,:), pointer, contiguous :: ostabl
        integer :: lunit, np6, ntemp, nnu, maxrec
        logical :: initialised=.false.
        real(dp) :: eosmet

    contains
        procedure :: init => init_osmet
        procedure :: add_line_opacities

    end type

contains
    
    !===============================================================================
    ! initialises osmet type.
    !===============================================================================
    SUBROUTINE init_osmet (self, absmet_file, absmet_big_end)
        class(osmet_type) :: self
        character(len=80), intent(in) :: absmet_file
        logical, intent(in), optional :: absmet_big_end
        integer :: lunit, iostat, maxp6, maxtem
        real(4) :: orig_ostabl(mp6, mtem)
        real(4), allocatable :: p6(:), temp(:)
        real(8) :: wnum
        integer :: nnu, nt, np
        !---------------------------------------------------------------------------
        print *, 'Reading absmet file: ' // trim(absmet_file)

        call os%check_file(absmet_file)
        if (absmet_big_end) then
            open(file=absmet_file, newunit=lunit, form='unformatted', convert='big_endian')
        else
            open(file=absmet_file, newunit=lunit, form='unformatted')
        endif

        read(lunit) !textatom
        read(lunit) !textatom
        read(lunit) !ispec
        read(lunit) !xiturb

        read(lunit) maxp6
        allocate(self%p6(maxp6), p6(maxp6))
        read(lunit) p6
        self%p6 = log10(dble(p6))
        
        read(lunit) maxtem
        allocate(self%temp(maxtem), temp(maxtem))
        read(lunit) temp
        self%temp = dble(temp)

        read(lunit)
        read(lunit)
        read(lunit)
        read(lunit)

        allocate(self%ostabl(maxp6, maxtem, mnu))
        self%ostabl = 0

        do nnu=1, mnu
            read(lunit, iostat=iostat) wnum, orig_ostabl
            if (iostat < 0) exit
            
            self%ostabl(:, :, mnu - nnu + 1) = orig_ostabl(:maxp6, :maxtem)
            self%nu(mnu - nnu + 1) = CC * dble(wnum)
        end do

        close(lunit)
        self%initialised = .true.
    END SUBROUTINE init_osmet

    !===============================================================================
    SUBROUTINE add_line_opacities(self, nu, p6, temp, ne, rho, abslte, sumsca)
        class(osmet_type), intent(inout) :: self
        real(8), dimension(:), intent(in) :: nu
        real(8), dimension(:,:,:), intent(in) :: p6, temp, ne, rho
        real(4), dimension(:,:,:,:), intent(inout) :: abslte, sumsca
        real(4), dimension(:,:,:), allocatable :: Tinv, sqrtT
        real(4) :: a1, b1, c1, d1, a2, b2, c2, d2
        real(4) :: ab1, cd1, ab2, cd2, abcd1, abcd2
        real(4) :: dx, dy, dz, ecnst, os
        real(dp) :: ecnst1, nu3, beta, epsln
        integer :: ix, iy, iz, i, j, k, f, nx, ny, nz, nf
        !---------------------------------------------------------------------------
        nx = size(abslte,1)
        ny = size(abslte,2)
        nz = size(abslte,3)
        nf = size(abslte,4)

        allocate(Tinv(nx,ny,nz))
        allocate(sqrtT(nx,ny,nz))
        Tinv = HH / (BK*temp)
        sqrtT = SQRT(temp)

        do f=1, nf
            iz = searchsorted(self%nu, nu(f))

            if (iz < 0) cycle
            
            dz = (nu(f) - self%nu(iz)) / (self%nu(iz+1) - self%nu(iz))
            nu3 = 1 / ( 20.6D0*(CC/nu(f))**3 )

            do k=1, nz
            do j=1, ny
            do i=1, nx
                ix = searchsorted(self%p6, p6(i,j,k))
                iy = searchsorted(self%temp, temp(i,j,k))

                if (ix<0 .or. iy<0) then
                    cycle
                endif

                a1 = self%ostabl(ix,   iy,   iz)
                b1 = self%ostabl(ix+1, iy,   iz)
                c1 = self%ostabl(ix+1, iy+1, iz)
                d1 = self%ostabl(ix,   iy+1, iz)
                a2 = self%ostabl(ix,   iy,   iz+1)
                b2 = self%ostabl(ix+1, iy,   iz+1)
                c2 = self%ostabl(ix+1, iy+1, iz+1)
                d2 = self%ostabl(ix,   iy+1, iz+1)
        
                dx = (p6  (i,j,k) - self%p6  (ix)) / (self%p6  (ix+1) - self%p6  (ix))
                dy = (temp(i,j,k) - self%temp(iy)) / (self%temp(iy+1) - self%temp(iy))
        
                ab1 = (b1 - a1) * dx + a1
                cd1 = (c1 - d1) * dx + d1
                ab2 = (b2 - a2) * dx + a2
                cd2 = (c2 - d2) * dx + d2
                
                abcd1 = exp(ln10*((cd1 - ab1) * dy + ab1))
                abcd2 = exp(ln10*((cd2 - ab2) * dy + ab2))
        
                ! last interpolation along frequencies is done on linear opacities
                os = (abcd2 - abcd1) * dz + abcd1
                os = os * rho(i,j,k)

                ecnst1 = sqrtT(i,j,k) * nu3
                beta = Tinv(i,j,k) * nu(f)
                ecnst = ecnst1*( self%eosmet / psemi(2,beta) + (1-self%eosmet) / psemi(1,beta) )
                epsln = ne(i,j,k) / (ecnst + ne(i,j,k))

                abslte(i,j,k,f) = abslte(i,j,k,f) + os*epsln
                sumsca(i,j,k,f) = sumsca(i,j,k,f) + os*(1 - epsln)
            end do
            end do
            end do
        end do
    END SUBROUTINE add_line_opacities

    !===============================================================================
    !  thermal p-FUNCTION for semic empirical collision rates
    !  reference: sobel'man - 'atomic physics ii.  '
    !===============================================================================
    real(8) FUNCTION psemi(z,b)
        integer, intent(in) :: z
        real(8), intent(in) :: b
        real(8), dimension(10), parameter :: betref = [-2.0,-1.699,-1.398,-1.0,-0.699,-0.398,0.0,0.301,0.602,1.0]
        real(8), dimension(10), parameter :: pnref  = [1.160,0.956,0.758,0.493,0.331,0.209,0.100,0.063,0.040,0.023]
        real(8), dimension(10), parameter :: pcref  = [1.160,0.977,0.788,0.554,0.403,0.290,0.214,0.201,0.200,0.200]
        real(8) :: bb
        integer :: i
        !---------------------------------------------------------------------------
        if (b .lt. 0.01) then
            !  limit of high temperature: kt >> e
            psemi = 0.27566 * (0.577 + log(b))
        elseif(b .gt. 0.01 .and. b .lt. 10.0)then
            !  intermediate temps (most important for eqm plasmas)
            !  spline interpolation onto logb grid:
            bb=log10(b)
            !if (z .eq. 1) psemi = splin(bb,betref,pnref,10,10)
            !if (z .gt. 1) psemi = splin(bb,betref,pcref,10,10)
            i = searchsorted(betref, bb)
            if (z .eq. 1) psemi = interp1d_linear(i,betref,pnref,bb)
            if (z .gt. 1) psemi = interp1d_linear(i,betref,pcref,bb)
        elseif(b .gt. 10.0)then
            !  limit of low temperature: kt << e
            if(z .eq. 1) psemi = 0.066 / sqrt(b)
            if(z .gt. 1) psemi = 0.200
        endif
        return
      end FUNCTION psemi
END MODULE osmet_mod