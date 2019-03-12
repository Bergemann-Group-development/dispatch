!===============================================================================
!===============================================================================
MODULE mpi_buffer_mod
  !USE trace_mod
  implicit none
  private
  type, public:: mpi_buffer_t
    integer:: n
    integer, allocatable:: buffer(:)
    integer(8):: position
  contains
    procedure:: init
    procedure:: reset
    procedure, private:: append_c
    procedure, private:: append_0
    procedure, private:: append_1
    procedure, private:: append_2
    procedure, private:: append_3
    generic, public:: append => append_c, append_1, append_2, append_3
    procedure, private:: read_c
    procedure, private:: read_0
    procedure, private:: read_1
    procedure, private:: read_2
    procedure, private:: read_3
    generic, public:: read => read_c, read_1, read_2, read_3
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, n, label)
  class (mpi_buffer_t):: self
  integer:: n
  character(len=*):: label
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%init')
  allocate (self%buffer(n))
  self%position = 1
  call self%append (label)
  !call trace%end
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE reset (self)
  class (mpi_buffer_t):: self
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%reset')
  self%position = 1
  !call trace%end
END SUBROUTINE reset

!===============================================================================
!===============================================================================
SUBROUTINE append_c (self, label)
  class (mpi_buffer_t):: self
  character(len=*):: label
  character, allocatable:: text(:)
  integer:: ni, nb
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  nb = len(label)
  ni = 1 + (nb-1)/4
  nb = ni*4
  allocate (text(nb))
  text = label
  self%buffer(self%position) = ni
  self%position = self%position + 1
  call mpi_buffer_copy (text, ni, self%buffer, self%position)
  self%position = self%position + ni
  deallocate (text)
  !call trace%end
END SUBROUTINE append_c

!===============================================================================
!===============================================================================
SUBROUTINE append_0 (self, a)
  class (mpi_buffer_t):: self
  real:: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  self%buffer(self%position) = 1
  self%position = self%position + 1
  call mpi_buffer_copy (a, 1, 1, self%buffer, self%position)
  self%position = self%position + 1
  !call trace%end
END SUBROUTINE append_0

!===============================================================================
!===============================================================================
SUBROUTINE append_1 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  self%buffer(self%position) = na
  self%position = self%position + 1
  call mpi_buffer_copy (a, 1, 1, self%buffer, self%position)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE append_1

!===============================================================================
!===============================================================================
SUBROUTINE append_2 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  self%buffer(self%position) = na
  self%position = self%position + 1
  call mpi_buffer_copy (a, 1, na, self%buffer, self%position)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE append_2

!===============================================================================
!===============================================================================
SUBROUTINE append_3 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%append_c')
  na = product(shape(a))
  self%buffer(self%position) = na
  self%position = self%position + 1
  call mpi_buffer_copy (a, 1, na, self%buffer, self%position)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE append_3

!===============================================================================
!===============================================================================
SUBROUTINE read_c (self, label)
  class (mpi_buffer_t):: self
  character(len=*):: label
  integer:: ni, nb
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  nb = len(label)
  ni = 1 + (nb-1)/4
  nb = ni*4
  !allocate (a(nb))
  self%position = self%position + 1
  !call mpi_buffer_copy (self%buffer, self%position, a, 1, ni)
  self%position = self%position + ni
  !deallocate (a)
  !call trace%end
END SUBROUTINE read_c

!===============================================================================
!===============================================================================
SUBROUTINE read_0 (self, a)
  class (mpi_buffer_t):: self
  real:: a
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  self%position = self%position + 1
  call mpi_buffer_copy (self%buffer, self%position, a, 1, 1)
  self%position = self%position + 1
  !call trace%end
END SUBROUTINE read_0

!===============================================================================
!===============================================================================
SUBROUTINE read_1 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  self%position = self%position + 1
  call mpi_buffer_copy (self%buffer, self%position, a, 1, na)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE read_1

!===============================================================================
!===============================================================================
SUBROUTINE read_2 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  self%buffer(self%position) = na
  self%position = self%position + 1
  call mpi_buffer_copy (self%buffer, self%position, a, 1, na)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE read_2

!===============================================================================
!===============================================================================
SUBROUTINE read_3 (self, a)
  class (mpi_buffer_t):: self
  real:: a(:,:,:)
  integer:: na
  !-----------------------------------------------------------------------------
  !call trace%begin ('mpi_buffer_t%read_c')
  na = product(shape(a))
  self%position = self%position + 1
  call mpi_buffer_copy (self%buffer, self%position, a, 1, na)
  self%position = self%position + na
  !call trace%end
END SUBROUTINE read_3

END MODULE mpi_buffer_mod

!===============================================================================
!> Copy n words from a to b, starting at position o
!===============================================================================
SUBROUTINE mpi_buffer_copy (a, n, b, o)
  integer:: a(:), n, b(:), o
  !-----------------------------------------------------------------------------
  b(o:o+n-1) = a(1:n)
END SUBROUTINE mpi_buffer_copy
