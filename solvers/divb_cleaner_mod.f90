!===============================================================================
!> Module to clear div(B) from the interiors of patches on a task list
!===============================================================================
MODULE divb_cleaner_mod
    USE io_mod
    USE trace_mod
    USE timer_mod
    USE kinds_mod
    USE task_mod
    USE patch_mod
    USE link_mod
    USE poisson_mod
    USE list_mod
    USE omp_mod
    implicit none
    private
    type, public:: divb_cleaner_t
    contains
      procedure:: init
      procedure:: update
    end type
    logical, save:: on=.false.
    integer, save:: verbose=0
    integer, save:: max_iter=5
    type(divb_cleaner_t), public:: divb_cleaner
  CONTAINS
  
  !===============================================================================
  !> Read namelist
  !===============================================================================
  SUBROUTINE init (self)
    class(divb_cleaner_t):: self
    integer:: iostat
    logical, save:: first_time=.true.
    namelist /divb_cleaner_params/ on, verbose, max_iter
    !-----------------------------------------------------------------------------
    if (first_time) then
      call io%lock%set ('divb_cleaner_t')
      rewind (io%input)
      read (io%input, divb_cleaner_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, divb_cleaner_params)
      end if
      write (io%namelist, divb_cleaner_params)
      call poisson%init ()
      first_time = .false.
      call io%lock%unset ('divb_cleaner_t')
    end if
  END SUBROUTINE init
  
  !===============================================================================
  !> Update the magnetic field in all patches in the task list, enforcing
  !> div(B)=0 internally -- guard zones are taken care of by the MHD solver
  !===============================================================================
  SUBROUTINE update (self, task_list)
    class(divb_cleaner_t):: self
    class(list_t):: task_list
    !.............................................................................
    class(link_t), pointer:: link
    class(task_t), pointer:: patch
    real(var_k), pointer, dimension(:,:,:,:):: bb
    real(var_k), pointer, dimension(:,:,:,:), contiguous:: dbb
    real(var_k), pointer, dimension(:,:,:), contiguous:: src, phi
    integer:: iter, i
    !-----------------------------------------------------------------------------
    if (.not. on) return
              if (verbose >= 0) then
                call io%header ('-cleaning div(B)', left=25)
                call timer%tic ()
              end if
    link => task_list%head
    i = 0
    do while (associated(link))
      !---------------------------------------------------------------------------
      ! Distribute work over OpenMP threads
      !---------------------------------------------------------------------------
      if (mod(i,omp%nthreads) == omp%thread) then
        call link%update ()
        patch => link%task
        select type (patch)
        class is (patch_t)
          bb => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%it,1)
          allocate (dbb(patch%gn(1),patch%gn(2),patch%gn(3),3))
          allocate (src(patch%gn(1),patch%gn(2),patch%gn(3)))
          allocate (phi(patch%gn(1),patch%gn(2),patch%gn(3)))
          do iter=1,max_iter
            call div (patch, bb, src)
            if (verbose > 0) &
              write (stdout,'(a,i6,i3,1p,g12.3)') &
                'divb_cleaner: id,iter,divb=', patch%id, iter, maxval(abs(src))
            call poisson%sor (patch, src, phi, src, patch%li-1, patch%ui+1)
            call grad (patch, phi, dbb)
            bb = bb - dbb
          end do
          deallocate (dbb, src, phi)
        end select
      end if
      link => link%next
    end do
              if (verbose >= 0) then
                call timer%toc ()
                call io%header ('-end cleaning div(B)', left=25)
              end if
  END SUBROUTINE update
  
  !===============================================================================
  !> Compute div(B) from up-staggered differences
  !===============================================================================
  SUBROUTINE div (patch, bb, divb)
    class(patch_t):: patch
    real(var_k), pointer, dimension(:,:,:,:):: bb
    real(var_k), pointer, dimension(:,:,:):: divb
    !.............................................................................
    real(8):: cs(3)
    integer:: l(3), u(3), i1, i2, i3
    !-----------------------------------------------------------------------------
    cs = 1d0/patch%ds
    l = patch%li
    u = patch%ui
    divb = 0.0
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      divb(i1,i2,i3) = (bb(i1+1,i2,i3,1) - bb(i1,i2,i3,1))*cs(1) &
                     + (bb(i1,i2+1,i3,2) - bb(i1,i2,i3,2))*cs(2) &
                     + (bb(i1,i2,i3+1,3) - bb(i1,i2,i3,3))*cs(3)
    end do
    end do
    end do
  END SUBROUTINE div
  
  !===============================================================================
  !> Compute grad(phi) from down-staggered differences
  !===============================================================================
  SUBROUTINE grad (patch, phi, grad_phi)
    class(patch_t):: patch
    real(var_k), pointer, dimension(:,:,:):: phi
    real(var_k), pointer, dimension(:,:,:,:):: grad_phi
    !.............................................................................
    real(8):: cs(3)
    integer:: l(3), u(3), i1, i2, i3
    !-----------------------------------------------------------------------------
    cs = 1d0/patch%ds
    l = patch%li
    u = patch%ui
    grad_phi = 0.0
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      grad_phi(i1,i2,i3,1) = (phi(i1,i2,i3) - phi(i1-1,i2,i3))*cs(1)
      grad_phi(i1,i2,i3,2) = (phi(i1,i2,i3) - phi(i1,i2-1,i3))*cs(2)
      grad_phi(i1,i2,i3,3) = (phi(i1,i2,i3) - phi(i1,i2,i3-1))*cs(3)
    end do
    end do
    end do
  END SUBROUTINE grad
  
  END MODULE divb_cleaner_mod
  