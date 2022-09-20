  !===============================================================================
!===============================================================================
MODULE mpi_buffer_mod
  !USE trace_mod
  implicit none
  private
  type, public:: mpi_buffer_t
    integer(8):: n=0
    integer(8):: position=0
    integer(8):: n_buffer=0
    logical:: dryrun=.false.
    integer, pointer, contiguous:: buffer(:)
  contains
    procedure:: init
    procedure:: reset
    procedure, private:: append_c
    procedure, private:: append_0
    procedure, private:: append_1
    procedure, private:: append_0d
    procedure, private:: append_0i
    procedure, private:: append_0id
    procedure, private:: append_1i
    procedure, private:: append_1id
    procedure, private:: append_1d
    procedure, private:: append_1q
    procedure, private:: append_2
    procedure, private:: append_2d
    procedure, private:: append_2i
    procedure, private:: append_2q
    procedure, private:: append_3
    procedure, private:: append_3i
    procedure, private:: append_3i2
    procedure, private:: append_4
    generic, public:: append => append_0, append_1, append_2, append_3, &
      append_c, append_1i, append_1q, append_1d, append_2q, append_0d, &
      append_0i, append_0id, append_1id, append_3i, append_3i2, append_4, &
      append_2i, append_2d
    procedure, private:: read_c
    procedure, private:: read_0
    procedure, private:: read_1
    procedure, private:: read_0d
    procedure, private:: read_0i
    procedure, private:: read_0id
    procedure, private:: read_1i
    procedure, private:: read_1id
    procedure, private:: read_1d
    procedure, private:: read_1q
    procedure, private:: read_2
    procedure, private:: read_2d
    procedure, private:: read_2q
    procedure, private:: read_3
    procedure, private:: read_3i
    procedure, private:: read_3i2
    procedure, private:: read_4
    generic, public:: read => read_0, read_1, read_2, read_3, &
      read_c, read_1i, read_1q, read_1d, read_2q, read_0d, &
      read_0i, read_0id, read_1id, read_3i, read_3i2, read_4, read_2d
  end type
CONTAINS

!===============================================================================
!> Initialize a buffer and optionally push a label onto it
!===============================================================================
SUBROUTINE init (self, n, label)
  class (mpi_buffer_t):: self
  integer(8):: n
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%init')
  allocate (self%buffer(n))
  self%n = n
  self%position = 1
  !call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate buffer and clear counter
!===============================================================================
SUBROUTINE reset (self)
  class (mpi_buffer_t):: self
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%reset')
  if (associated(self%buffer)) deallocate(self%buffer)
  self%position = 0
  !call trace%end()
END SUBROUTINE reset

!===============================================================================
!> Push a text label onto the buffer (must know the size when retrieving)
!===============================================================================
SUBROUTINE append_c (self, label)
  class (mpi_buffer_t):: self
  character(len=*):: label
  character, allocatable:: a(:)
  integer:: na, nb
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  nb = len(label)
  na = 1 + (nb-1)/4
  nb = na*4
  if (self%dryrun) then
    self%n = self%n + na
  else
    allocate (a(nb))
    a = label
    call to_buffer (a, na, self%buffer, self%n, self%position)
    deallocate (a)
  end if
  !call trace%end()
END SUBROUTINE append_c

!===============================================================================
!===============================================================================
SUBROUTINE append_0 (self, a)
  class (mpi_buffer_t):: self
  real:: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  if (self%dryrun) then
    self%n = self%n + 1
  else
    call to_buffer1 (a, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_0

!===============================================================================
!===============================================================================
SUBROUTINE append_0i (self, a)
  class (mpi_buffer_t):: self
  integer:: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  if (self%dryrun) then
    self%n = self%n + 1
  else
    call to_buffer1 (a, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_0i

!===============================================================================
!===============================================================================
SUBROUTINE append_0id (self, a)
  class (mpi_buffer_t):: self
  integer(8):: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  if (self%dryrun) then
    self%n = self%n + 2
  else
    call to_buffer (a, 2, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_0id

!===============================================================================
!===============================================================================
SUBROUTINE append_0d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  if (self%dryrun) then
    self%n = self%n + 2
  else
    call to_buffer (a, 2, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_0d

!===============================================================================
!===============================================================================
SUBROUTINE append_1 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = size(a)
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_1

!===============================================================================
!===============================================================================
SUBROUTINE append_1d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = size(a)*2
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_1d

!===============================================================================
!===============================================================================
SUBROUTINE append_1q (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = size(a)/2
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_1q

!===============================================================================
!===============================================================================
SUBROUTINE append_1i (self, a)
  class (mpi_buffer_t):: self
  integer:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = size(a)
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_1i


!===============================================================================
!===============================================================================
SUBROUTINE append_1id (self, a)
  class (mpi_buffer_t):: self
  integer(8):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = size(a)*2
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_1id

!===============================================================================
!===============================================================================
SUBROUTINE append_2 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_2

!===============================================================================
!===============================================================================
SUBROUTINE append_2d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  na = product(shape(a))*2
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
END SUBROUTINE append_2d

!===============================================================================
!===============================================================================
SUBROUTINE append_2i (self, a)
  class (mpi_buffer_t):: self
  integer:: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_2i

!===============================================================================
!> Convert integer(2) to integer, to avoid Python problems
!===============================================================================
SUBROUTINE append_2q (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_2q

!===============================================================================
!===============================================================================
SUBROUTINE append_3 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_3

!===============================================================================
!===============================================================================
SUBROUTINE append_3i (self, a)
  class (mpi_buffer_t):: self
  integer:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_3i

!===============================================================================
!===============================================================================
SUBROUTINE append_3i2 (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_3i2

!===============================================================================
!===============================================================================
SUBROUTINE append_4 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  if (self%dryrun) then
    self%n = self%n + na
  else
    call to_buffer (a, na, self%buffer, self%n, self%position)
  end if
  !call trace%end()
END SUBROUTINE append_4

!===============================================================================
!===============================================================================
SUBROUTINE read_c (self, a)
  class (mpi_buffer_t):: self
  character(len=*):: a
  integer:: na, nb
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  nb = len(a)
  na = 1 + (nb-1)/4
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_c

!===============================================================================
!===============================================================================
SUBROUTINE read_0 (self, a)
  class (mpi_buffer_t):: self
  real:: a
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  call from_buffer (a, 1, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_0

!===============================================================================
!===============================================================================
SUBROUTINE read_0i (self, a)
  class (mpi_buffer_t):: self
  integer:: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  call from_buffer1 (a, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_0i

!===============================================================================
!===============================================================================
SUBROUTINE read_0id (self, a)
  class (mpi_buffer_t):: self
  integer(8):: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  call from_buffer (a, 2, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_0id

!===============================================================================
!===============================================================================
SUBROUTINE read_0d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  call from_buffer (a, 2, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_0d

!===============================================================================
!===============================================================================
SUBROUTINE read_1 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = size(a)
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_1

!===============================================================================
!===============================================================================
SUBROUTINE read_1d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = size(a)*2
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_1d

!===============================================================================
!===============================================================================
SUBROUTINE read_1q (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = size(a)/2
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_1q

!===============================================================================
!===============================================================================
SUBROUTINE read_1i (self, a)
  class (mpi_buffer_t):: self
  integer:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = size(a)
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_1i

!===============================================================================
!===============================================================================
SUBROUTINE read_1id (self, a)
  class (mpi_buffer_t):: self
  integer(8):: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = size(a)*2
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_1id

!===============================================================================
!===============================================================================
SUBROUTINE read_2 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_2

!===============================================================================
!===============================================================================
SUBROUTINE read_2d (self, a)
  class (mpi_buffer_t):: self
  real(8):: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  na = product(shape(a))*2
  call from_buffer (a, na, self%buffer, self%n, self%position)
END SUBROUTINE read_2d

!===============================================================================
!===============================================================================
SUBROUTINE read_2q (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))/2
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_2q

!===============================================================================
!===============================================================================
SUBROUTINE read_3 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_3

!===============================================================================
!===============================================================================
SUBROUTINE read_3i (self, a)
  class (mpi_buffer_t):: self
  integer:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_3i

!===============================================================================
!===============================================================================
SUBROUTINE read_3i2 (self, a)
  class (mpi_buffer_t):: self
  integer(2):: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_3i2

!===============================================================================
!===============================================================================
SUBROUTINE read_4 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  call from_buffer (a, na, self%buffer, self%n, self%position)
  !call trace%end()
END SUBROUTINE read_4

END MODULE mpi_buffer_mod

!===============================================================================
!> Copy n words from a to b, starting at position o
!===============================================================================
SUBROUTINE to_buffer (a, n, b, m, o)
  integer:: a(n), n, m, b(m), o
  !-----------------------------------------------------------------------------
  if (n > 0) then
    b(o:o+n-1) = a(1:n)
    o = o+n
  end if
END SUBROUTINE to_buffer

SUBROUTINE to_buffer1 (a, b, m, o)
  integer:: a, m, b(m), o
  !-----------------------------------------------------------------------------
  b(o)=a
  o = o+1
END SUBROUTINE to_buffer1

SUBROUTINE from_buffer (a, n, b, m, o)
  integer:: a(n), n, b(m), m, o
  !-----------------------------------------------------------------------------
  a(1:n) = b(o:o+n-1)
  o = o+n
END SUBROUTINE from_buffer

SUBROUTINE from_buffer1 (a, n, b, o)
  integer:: a, n, b(n), o
  !-----------------------------------------------------------------------------
  a = b(o)
  o = o+1
END SUBROUTINE from_buffer1
