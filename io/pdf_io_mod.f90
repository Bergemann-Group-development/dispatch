!===============================================================================
!> Module to periodically output a density probability density function (PDF).
!> The procedure assumes that 'out_time' is larger than the largest time step.
!===============================================================================
MODULE pdf_io_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE link_mod
  USE task_mod
  USE counters_mod
  USE mpi_mod
  USE bits_mod
  implicit none
  private
  !.............................................................................
  integer, parameter:: n_rms=6
  type, public:: pdf_io_t
    logical:: on=.false.
    integer:: iout=0
    real(8):: time, out_time
    real:: rms(n_rms)=0.0
    type(counters_t):: counters
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  !.............................................................................
  real, dimension(:), allocatable:: bins, counts
  real(8):: out_time=0.1_8
  real:: bin_min=-12.0, bin_max=10.0, bin_delta
  integer:: nbins=111
  integer:: verbose=0
  logical:: on=.false.
  type(pdf_io_t), public:: pdf_io
CONTAINS

!===============================================================================
!> Initialize parameters for PDF output
!===============================================================================
SUBROUTINE init (self, patch)
  class(pdf_io_t):: self
  class(patch_t):: patch
  !.............................................................................
  integer:: i, iostat
  logical:: first_time=.true.
  namelist /pdf_io_params/ on, verbose, out_time, nbins, bin_min, bin_max
  !-----------------------------------------------------------------------------
  call trace%begin ('pdf_io_t%init')
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      rewind (io%input)
      read (io%input, pdf_io_params, iostat=iostat)
      write (io%output, pdf_io_params)
      if (on) then
        allocate (bins(nbins), counts(nbins))
        bin_delta = (bin_max-bin_min)/(nbins-1.0)
        do i=1,nbins
          bins(i) = bin_min + (i-1)*bin_delta
          counts(i) = 0.0
        end do
      end if
      first_time = .false.
    end if
    !$omp end critical (input_cr)
  end if
  self%on = on
  if (on) then
    self%out_time = out_time
    if (patch%pdf_next==0d0) &
      patch%pdf_next = out_time
    call self%counters%init
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> PDF update, to be called when a patch first passes out_next
!===============================================================================
SUBROUTINE update (self, patch)
  class(pdf_io_t):: self
  class(patch_t):: patch
  !.............................................................................
  character(len=80):: filename
  integer:: ix, iy, iz, l(3), u(3), i, ibin, jt(2), iout, count
  real:: dv, bin, pt(2), d, px, py, pz, ux, uy, uz, u2, rms(n_rms)
  real(8):: v_tot, d_tot, px_tot, py_tot, pz_tot, u2_tot
  logical:: deallocated
  !-----------------------------------------------------------------------------
  if (.not.on) return
  if (patch%is_set(bits%remove)) then
    write (stdout,*) 'WARNING: tried to count PDF of removed', patch%id
    return
  end if
  if (patch%is_set(bits%virtual)) then
    write (stderr,*) 'ERROR: tried to count PDF of virtual', patch%id
    return
  end if
  call trace%begin ('pdf_io_t%update')
  if (patch%time > patch%pdf_next) then
    call patch%time_interval (patch%pdf_next, jt, pt)
    deallocated = .not.allocated(patch%irefine)
    if (deallocated) &
      allocate (patch%irefine(patch%gn(1),patch%gn(2),patch%gn(3)))
    call make_irefine (patch)
    dv = product(patch%mesh%d)
    d_tot  = 0.0
    px_tot = 0.0
    py_tot = 0.0
    pz_tot = 0.0
    u2_tot = 0.0
    v_tot  = 0.0
    l = patch%mesh%li
    u = patch%mesh%ui
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      if (patch%irefine(ix,iy,iz) == 0) then
        d =  patch%mem(ix,iy,iz,patch%idx%d ,jt(1),1)*pt(1) &
           + patch%mem(ix,iy,iz,patch%idx%d ,jt(2),1)*pt(2)
        px = patch%mem(ix,iy,iz,patch%idx%px,jt(1),1)*pt(1) &
           + patch%mem(ix,iy,iz,patch%idx%px,jt(2),1)*pt(2)
        py = patch%mem(ix,iy,iz,patch%idx%py,jt(1),1)*pt(1) &
           + patch%mem(ix,iy,iz,patch%idx%py,jt(2),1)*pt(2)
        pz = patch%mem(ix,iy,iz,patch%idx%pz,jt(1),1)*pt(1) &
           + patch%mem(ix,iy,iz,patch%idx%pz,jt(2),1)*pt(2)
        ux = px/d
        uy = py/d
        uz = pz/d
        u2 = ux**2 + uy**2 + uz**2
        px_tot = px_tot + px*dv
        py_tot = py_tot + py*dv
        pz_tot = pz_tot + pz*dv
        u2_tot = u2_tot + u2*dv
        d_tot  = d_tot  +  d*dv
        v_tot  = v_tot  +    dv
        bin = (alog(d) - bin_min)/bin_delta
        ibin = 1.0+bin
        ibin = max(1,min(ibin,nbins))
        counts(ibin) = counts(ibin) + dv
      end if
    end do
    end do
    end do
    rms = [v_tot, d_tot, px_tot, py_tot, pz_tot, u2_tot]
    do i=1,n_rms
      !$omp atomic update
      self%rms(i) = self%rms(i) + rms(i)
    end do
    if (deallocated) &
      deallocate (patch%irefine)
    !---------------------------------------------------------------------------
    ! Decrement the count of tasks still remaining to be accumulated.  iout is
    ! the output index (0 at time=0.0) and the counter index, which needs to be
    ! non-zero, is iout+1.  The first task that passes a give pdf_next refers to
    ! a non-existing counter, which causes the counter to be created, with initial
    ! count = io%nwrite-1.  
    !---------------------------------------------------------------------------
    !$omp critical (update_counters_cr)
    iout = nint(patch%pdf_next/pdf_io%out_time)
    call self%counters%update (iout+1, io%nwrite, -1, count)
    if (verbose > 0) &
      write(stdout,'(a,i6,3f12.6,4i4,i5,2x,a)') &
        'pdf_io_t%update: id, time, pdf_next =', patch%id, patch%time, &
        patch%time-patch%dtime, patch%pdf_next, count, self%iout, iout, &
        patch%level, patch%n_nbors, 'output'
    !---------------------------------------------------------------------------
    ! Call the output procedure when all tasks have been accumulated.  All tasks
    ! have their pdf_next incremented after their data are accumulated.
    !---------------------------------------------------------------------------
    if (count==0) then
      self%time = patch%pdf_next
      call self%output
    end if
    patch%pdf_next = (iout+1)*out_time
    !$omp end critical (update_counters_cr)
  end if
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> PDF output, delayed until all patches have been processed.  The PDF snapshot
!> number is kept as pdf_io%iout, and since we only use one counts(:) array, we
!> need to assume that this is the same for all patches.
!===============================================================================
SUBROUTINE output (self)
  class(pdf_io_t):: self
  !.............................................................................
  character(len=80):: file
  integer:: pdf_iout
  integer, parameter:: ioformat=1
  !-----------------------------------------------------------------------------
  call trace%begin ('pdf_io_t%output')
  !---------------------------------------------------------------------------
  ! Increment the output counter, keeping track of the no of snapshots written.
  ! The counter corresponding to the pending snapshot is (or has been) created
  ! when the first task sets pdf_next to that output time.
  !---------------------------------------------------------------------------
  self%iout = self%iout+1
  write (file,'(a,"pdf_",i5.5,"_",i5.5,".dat")') &
    trim(io%outputname), self%iout, mpi%rank
  open  (io_unit%tmp, file=file, form='unformatted', status='unknown')
  write (io_unit%tmp) ioformat, nbins
  write (io_unit%tmp) self%time
  write (io_unit%tmp) bins
  write (io_unit%tmp) counts
  write (io_unit%tmp) self%rms
  close (io_unit%tmp)
  if (verbose > -1) &
    write (stdout,*) trim(file)
  counts(:) = 0.0
  self%rms(:) = 0.0
  call trace%end()
END SUBROUTINE output

!===============================================================================
!> Make a map of where the patch is refined
!===============================================================================
SUBROUTINE make_irefine (patch)
  class(patch_t):: patch
  class(link_t), pointer:: nbor
  class(task_t), pointer:: nbpatch
  integer:: l(3), u(3)
  real(8):: dist(3)
  !-----------------------------------------------------------------------------
  patch%irefine = 0
  nbor => patch%link%nbor
  do while (associated(nbor))
    nbpatch => nbor%task
    if (nbpatch%level > patch%level+1) then
      write (stderr,*) 'WARNING: nbpatch w/o support =', nbpatch%id, patch%id
    else if (nbpatch%level > patch%level) then
      select type (nbpatch)
      class is (patch_t)
        if (patch%contains(nbpatch)) then
          dist = (nbpatch%size-nbpatch%ds)/2d0
          l = patch%index_only (nbpatch%position-dist)
          u = patch%index_only (nbpatch%position+dist)
          if (verbose > 1) &
            write(stdout,*) 'irefine: nbpatch =', &
              nbpatch%id, nbpatch%level, l, u, product(u-l+1)
          patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 1
        end if
      end select
    end if
    nbor => nbor%next
  end do
  if (verbose > 0) then
    l = patch%li
    u = patch%ui
    write(stdout,*) 'make_irefine: unrefined =', &
      sum(1-patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3))), product(u-l+1)
  end if
END SUBROUTINE make_irefine

END MODULE pdf_io_mod
