!===============================================================================
!> Flux currection module
!===============================================================================
MODULE flux_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  implicit none
  private
  type, public:: flux_t
  contains
    procedure:: init
    procedure:: corrections
    procedure:: correct
  end type
  integer, parameter:: nf=3
  logical, save:: on=.false.
  integer, save:: verbose=0
  type(flux_t), public:: flux
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self)
  class(flux_t):: self
  integer:: iostat
  namelist /flux_params/ on, verbose
  logical, save:: first_time=.true.
  !.............................................................................
  call trace%begin ('flux_t%init')
  if (first_time) then
    call io%lock%set
    if (first_time) then
      rewind (io_unit%input)
      read (io%input, flux_params, iostat=iostat)
      call io%namelist_warning ('flux_params', iostat)
      write (io%namelist, flux_params)
       first_time = .false.
    end if
    call io%lock%unset
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Select nbors that share a cube face with the give patch
!===============================================================================
SUBROUTINE corrections (self, patch)
  class(flux_t):: self
  class(patch_t):: patch
  !.............................................................................
  class(nbor_t), pointer:: nbor
  real(8):: rr(3), eps=1e-2
  integer:: ii(3), dir, face
  !-----------------------------------------------------------------------------
  call trace%begin ('flux_t%corrections')
  call self%init
  nbor => patch%link%nbor
  do while (associated(nbor))
    rr = patch%distance(nbor%task)/(0.5d0*(patch%size+nbor%task%size))
    if (verbose > 1) &
      write(io%output,'(a,2i6,2x,3f8.3)') &
        'flux_t%correct: patch, nbor, rr =', patch%id, nbor%task%id, rr
    !---------------------------------------------------------------------------
    ! Only patches that are outside in exactly one direction have overlap with
    ! a cube face
    !---------------------------------------------------------------------------
    if (count(abs(abs(rr)-1d0) < eps) == 1) then
      do dir=1,3
        do face=-1,1,2
          if (abs(rr(dir)-face) < eps) then
            associate (nbpatch => nbor%task)
            select type (nbpatch)
            class is (patch_t)
              call self%correct (patch, nbpatch, dir, face)
            end select
            end associate
          end if
        end do
      end do
    end if
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE corrections

!===============================================================================
!> Correct the face values, based on the flux values
!===============================================================================
SUBROUTINE correct (self, patch, nbpatch, idir, iface)
  class(flux_t):: self
  class(patch_t):: patch, nbpatch
  character(len=2):: dir(3) = ['x-', 'y-', 'z-']
  character(len=5):: face(-1:1) = ['left ', '     ', 'right']
  integer:: idir, iface
  !.............................................................................
  class(link_t), pointer:: nbor
  real, dimension(:,:,:,:,:), allocatable:: mem
  !-----------------------------------------------------------------------------
  ! Allocate working memory
  !-----------------------------------------------------------------------------
  call trace%begin ('flux_t%correct')
  if (verbose > 0) &
    write(stdout,'(a,2i6,3x,a)') &
      'flux_t%correct: ', patch%id, nbpatch%id, dir(idir)//face(iface)
  select case (idir)
  case (1)
    allocate (mem(patch%n(2),patch%n(3),nf,patch%nv,2))
  case (2)
    allocate (mem(patch%n(1),patch%n(3),nf,patch%nv,2))
  case (3)
    allocate (mem(patch%n(1),patch%n(2),nf,patch%nv,2))
  end select
  !-----------------------------------------------------------------------------
  ! Gather face values into mem, update, and scatter back
  !-----------------------------------------------------------------------------
  if (on) then
    call gather_scatter (self, patch, nbpatch, idir, iface, mem, .true.)
    call update (mem)
    call gather_scatter (self, patch, nbpatch, idir, iface, mem, .false.)
  end if
  !-----------------------------------------------------------------------------
  deallocate (mem)
  call trace%end ()
END SUBROUTINE correct

!===============================================================================
!> Gather or scatter face values to / from working mem
!===============================================================================
SUBROUTINE gather_scatter (self, patch, nbpatch, dir, face, mem, gather)
  class(flux_t):: self
  class(patch_t):: patch, nbpatch
  integer:: dir, face
  real:: mem(:,:,:,:,:)
  logical:: gather
  !.............................................................................
  integer:: i, i0, j0, j1, i1, i2, i3, iv, iw
  !-----------------------------------------------------------------------------
  ! Cast nbpatch to patch_t
  !-----------------------------------------------------------------------------
  call trace%begin ('flux_t%gather')
  !-----------------------------------------------------------------------------
  ! Start index and index direction
  !-----------------------------------------------------------------------------
  select case (face)
  case (+1)
    i0 = patch%mesh(dir)%ui-1
    j1 = +1
  case (-1)
    i0 = patch%mesh(dir)%li+1
    j1 = -1
  end select
  !-----------------------------------------------------------------------------
  ! Pick up values for the chosen face, from inside to outside
  !-----------------------------------------------------------------------------
  select case (dir)
  case (1)
    do iv=1,patch%nv
      do i3=patch%mesh(3)%li,patch%mesh(3)%ui
        i1 = i0
        do i=1,nf
          if (gather) then
            do iw=1,2
            do i2=patch%mesh(2)%li,patch%mesh(2)%ui
              mem(i2,i3,i,iv,iw) = patch%mem(i1,i2,i3,iv,patch%it,iw)
            end do
            end do
          else
            do i2=patch%mesh(2)%li,patch%mesh(2)%ui
              patch%mem(i1,i2,i3,iv,patch%it,1) = mem(i2,i3,i,iv,i)
            end do
          end if
          i1 = i1+j1
        end do
      end do
    end do
  case (2)
    do iv=1,patch%nv
      do i3=patch%mesh(3)%li,patch%mesh(3)%ui
        i2 = i0
        do i=1,nf
          if (gather) then
            do iw=1,2
              do i1=patch%mesh(1)%li,patch%mesh(1)%ui
                mem(i1,i3,i,iv,iw) = patch%mem(i1,i2,i3,iv,patch%it,iw)
              end do
            end do
          else
            do i1=patch%mesh(1)%li,patch%mesh(1)%ui
              patch%mem(i1,i2,i3,iv,patch%it,1) = mem(i1,i3,i,iv,1)
            end do
          end if
          i2 = i2+j1
        end do
      end do
    end do
  case (3)
    do iv=1,patch%nv
      do i2=patch%mesh(2)%li,patch%mesh(2)%ui
        i3 = i0
        do i=1,nf
          if (gather) then
            do iw=1,2
              do i1=patch%mesh(1)%li,patch%mesh(1)%ui
                mem(i1,i2,i,iv,iw) = patch%mem(i1,i2,i3,iv,patch%it,iw)
              end do
            end do
          else
            do i1=patch%mesh(1)%li,patch%mesh(1)%ui
              patch%mem(i1,i2,i3,iv,patch%it,1) = mem(i1,i2,i,iv,1)
            end do
          end if
          i3 = i3+j1
        end do
      end do
    end do
  end select
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE gather_scatter

SUBROUTINE update (mem)
  real:: mem(:,:,:,:,:)
  !-----------------------------------------------------------------------------
END SUBROUTINE update

END MODULE flux_mod
