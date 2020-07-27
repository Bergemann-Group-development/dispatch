!===============================================================================
!> Module to produce a file format that allows fast reading with Python
!===============================================================================
MODULE particle_io_mod
  USE io_mod
  USE particles_mod
  implicit none
  private
  type, public:: particle_io_t
    class(gpatch_t), pointer:: patch
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  integer:: n = 0
  integer(8):: n_tot = 0
  integer(2), pointer, contiguous :: q(:,:)
  integer(8), pointer, contiguous :: id(:)
  real(8),    pointer, contiguous :: t(:)
  real,       pointer, contiguous :: r(:,:), v(:,:), vg(:,:), w(:), d(:)
  !$omp threadprivate (n,q,id,t,r,v,vg,w,d)
  integer:: m = 1000
  integer:: io_max = 0, id_max = 0
  integer, parameter:: EOF = -1, particle_io_key = 2
CONTAINS

!===============================================================================
!> Initialize module
!===============================================================================
SUBROUTINE output (self, patch, lb, ub)
  class(particle_io_t) :: self
  class(particles_t) :: particles
  class(gpatch_t), pointer:: patch
  integer:: n
  !-----------------------------------------------------------------------------
  n = 0
  do i4=1,n_bin
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    n = n + self%particles(i1,i2,i3,i4)%n
  end do
  end do
  end do
  end do
  allocate (q(n,4), r(n,3), v(n,3), vg(n,3), d(n), w(n), t(n), id(n))
  class(particle_io_t)   :: self
  type(particles_t)          :: pp
  !.............................................................................
  class(particle_t), pointer :: p
  integer                    :: j
  !-----------------------------------------------------------------------------
  p => pp%head%next
  do j=1,pp%n
    select type (p)
    type is (vparticle_t)
      n = n + 1
      q (n,:) = p%q
      r (n,:) = p%r
      v (n,:) = p%v
      vg(n,:) = p%vg
      w (n)   = p%w
      d (n)   = p%d
      t (n)   = p%t
      id(n)   = p%id
      if (n == m) then
        call self%output
      end if
    end select
    p => p%next
  end do
END SUBROUTINE update

!===============================================================================
!> Module to produce a file format that allows fast reading with Python
!===============================================================================
SUBROUTINE output (self, iout)
  class(particle_io_t):: self
  integer, optional:: iout
  !-----------------------------------------------------------------------------
  ! Write out info on dimensions -- before data to ensure it is up-to-date
  ! during live runs
  !-----------------------------------------------------------------------------
  if (present(iout)) then
    id_max = max(id_max,maxval(id))
    io_max = max(io_max,iout)
    rewind (io_unit%particle_i)
    write (io_unit%particle_i) id_max, io_max
    flush (io_unit%particle_i)
  end if
  !-----------------------------------------------------------------------------
  ! Write out buffered data
  !-----------------------------------------------------------------------------
  write (io_unit%particle_d) particle_io_key
  write (io_unit%particle_d) m
  write (io_unit%particle_d) q, r, v, vg, w, d, t, id
  write (io_unit%particle_d) EOF
  backspace (io_unit%particle_d)
  n_tot = n_tot + m
  if (mod(n_tot,1000*m) == 0) then
    write(stdout,'(a,f8.2,a)') 'particle_io:', n_tot*1e-6, ' Mp'
  end if
  n = 0
END SUBROUTINE output

END MODULE particle_io_mod
