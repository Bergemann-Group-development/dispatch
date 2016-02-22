!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  implicit none
  public
  !
  integer, private:: verbose=0
  type, private:: node_t
    type(scalar_t), pointer:: scalar=>null()
    type(node_t), pointer:: next=>null()
  end type
  type(node_t), private, pointer:: head=>null(), tail=>null()
  integer, private, save:: id=0
  !
  type, public:: scalar_t
    real, dimension(:,:,:), pointer:: f=>null()
    logical:: allocated=.false., temporary=.false.
    integer:: n(3), l(3), u(3)
    integer:: id
    real:: ds(3)
    character(len=32):: name = ''
  contains
    procedure:: init => init_scalar
    procedure:: allocate => allocate_scalar
    procedure:: tempalloc
    procedure:: tempdealloc
  end type
  integer, private:: nalloc=0
  interface assignment(=)
    module procedure scalar_scalar_assign, scalar_assign, scalar_assign_real, scalar_assign_int
  end interface
  interface operator(+)
    module procedure scalar_scalar_add, scalar_real_add, scalar_int_add
  end interface
  interface operator(**)
    module procedure scalar_power_int, scalar_power_real
  end interface
  interface operator(*)
    module procedure scalar_real_mul, real_scalar_mul, real8_scalar_mul, scalar_int_mul, scalar_scalar_mul
  end interface
  interface operator(/)
    module procedure scalar_scalar_divide
  end interface
  interface operator(-)
    module procedure scalar_scalar_sub, scalar_real_sub, scalar_int_sub, scalar_minus
  end interface
  interface exp
    module procedure scalar_exp
  end interface
  interface log
    module procedure scalar_log
  end interface
  interface sqrt
    module procedure scalar_sqrt
  end interface
  interface minval
    module procedure scalar_minval
  end interface
  interface maxval
    module procedure scalar_maxval
  end interface
  interface max
    module procedure scalar_max
  end interface
  public init_scalar
CONTAINS
!> =============================================================================
SUBROUTINE init_scalar (self, dims, ds, name)
  class(scalar_t):: self
  integer:: dims(3)
  real:: ds(3)
  character(len=*), optional:: name
  !.............................................................................
  self%n = dims
  self%l = 1
  self%u = dims
  self%ds = ds
  !$omp atomic
  id = id+1
  self%id = id
  self%name=''; if (present(name)) self%name=name
  call self%allocate
  self%allocated = .true.
END SUBROUTINE init_scalar
!> =============================================================================
SUBROUTINE scalar_scalar_assign (out, in)
  type(scalar_t), intent(inout):: out
  type(scalar_t), intent(in):: in
  if (verbose>0) print*,out%id,out%name,out%allocated
  call out%tempalloc (in,'scalar =')
  out%ds = in%ds
  out%n  = in%n
  out%l  = in%l
  out%u  = in%u
  out%f  = in%f
  call in%tempdealloc ('scalar =')
END SUBROUTINE scalar_scalar_assign
!> =============================================================================
SUBROUTINE scalar_assign_real (out, in)
  type(scalar_t), intent(inout):: out
  real, intent(in):: in
  out%f = in
END SUBROUTINE scalar_assign_real
!> =============================================================================
SUBROUTINE scalar_assign_int (out, in)
  type(scalar_t), intent(inout):: out
  integer, intent(in):: in
  out%f = in
END SUBROUTINE scalar_assign_int
!> =============================================================================
SUBROUTINE scalar_assign (out, in)
  type(scalar_t), intent(inout):: out
  real, dimension(:,:,:), intent(in), pointer:: in
  out%f = in
END SUBROUTINE scalar_assign
!> =============================================================================
FUNCTION scalar_minus (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%tempalloc(in)
  out%f = - in%f
  call in%tempdealloc
END FUNCTION scalar_minus
!> =============================================================================
FUNCTION scalar_maxval (in) RESULT (out)
  real:: out
  type(scalar_t), intent(in):: in
  out = maxval(in%f)
  call in%tempdealloc
END FUNCTION scalar_maxval
!> =============================================================================
FUNCTION scalar_max (in1, in2) RESULT (out)
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  type(scalar_t):: out
  call out%tempalloc (in1)
  out%f = max(in1%f,in2)
  call in1%tempdealloc
END FUNCTION scalar_max
!> =============================================================================
FUNCTION scalar_minval (in) RESULT (out)
  real:: out
  type(scalar_t), intent(in):: in
  out = minval(in%f)
  call in%tempdealloc
END FUNCTION scalar_minval
!> =============================================================================
FUNCTION scalar_exp (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%tempalloc(in)
  out%f = exp(in%f)
  call in%tempdealloc
END FUNCTION scalar_exp
!> =============================================================================
FUNCTION scalar_log (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%tempalloc(in)
  out%f = log(in%f)
  call in%tempdealloc
END FUNCTION scalar_log
!> =============================================================================
FUNCTION scalar_sqrt (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%tempalloc(in)
  out%f = sqrt(in%f)
  call in%tempdealloc
END FUNCTION scalar_sqrt

!> =============================================================================
FUNCTION scalar_real_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f + in2
  call in1%tempdealloc
END FUNCTION scalar_real_add
!> =============================================================================
FUNCTION scalar_int_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f + in2
  call in1%tempdealloc
END FUNCTION scalar_int_add
!> =============================================================================
FUNCTION scalar_scalar_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%tempalloc(in1)
  out%f = in1%f + in2%f
  call in1%tempdealloc ('scalar add')
  call in2%tempdealloc ('scalar add')
END FUNCTION scalar_scalar_add

!> =============================================================================
FUNCTION real8_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real(8), intent(in):: in1
  type(scalar_t), intent(in):: in2
  call out%tempalloc(in2)
  out%f = in1 * in2%f
  call in2%tempdealloc
END FUNCTION real8_scalar_mul
!> =============================================================================
FUNCTION real_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in1
  type(scalar_t), intent(in):: in2
  call out%tempalloc(in2)
  out%f = in1 * in2%f
  call in2%tempdealloc
END FUNCTION real_scalar_mul
!> =============================================================================
FUNCTION scalar_real_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in2
  type(scalar_t), intent(in):: in1
  call out%tempalloc(in1)
  out%f = in2 * in1%f
  call in1%tempdealloc
END FUNCTION scalar_real_mul
!> =============================================================================
FUNCTION scalar_int_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f * in2
  call in1%tempdealloc
END FUNCTION scalar_int_mul
!> =============================================================================
FUNCTION scalar_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%tempalloc(in1)
  out%f = in1%f * in2%f
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION scalar_scalar_mul
!> =============================================================================
FUNCTION scalar_scalar_divide (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%tempalloc(in1)
  out%f = in1%f / in2%f
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION scalar_scalar_divide

!> =============================================================================
FUNCTION scalar_power_int (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f**in2
  call in1%tempdealloc
END FUNCTION scalar_power_int
!> =============================================================================
FUNCTION scalar_power_real (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f**in2
  call in1%tempdealloc
END FUNCTION scalar_power_real

!> =============================================================================
FUNCTION scalar_real_sub (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f - in2
  call in1%tempdealloc
END FUNCTION scalar_real_sub
!> =============================================================================
FUNCTION scalar_int_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%tempalloc(in1)
  out%f = in1%f - in2
  call in1%tempdealloc
END FUNCTION scalar_int_sub
!> =============================================================================
FUNCTION scalar_scalar_sub (in1, in2) RESULT (out)
  type(scalar_t), intent(in):: in1, in2
  type(scalar_t):: out
  call out%tempalloc(in1,'scalar sub')
  out%f = in1%f - in2%f
  call in1%tempdealloc ('scalar sub')
  call in2%tempdealloc ('scalar sub')
END FUNCTION scalar_scalar_sub

!> =============================================================================
SUBROUTINE allocate_scalar (self)
  class(scalar_t):: self
  !.............................................................................
  if (self%allocated) return
  !$omp atomic
  id = id+1
  self%id = id
  allocate (self%f(self%n(1), self%n(2), self%n(3)))
  self%allocated = .true.
END SUBROUTINE allocate_scalar

!> =============================================================================
SUBROUTINE tempalloc (self, in, label)
  class(scalar_t), target:: self
  type(scalar_t):: in
  character(len=*), optional:: label
  type(node_t), pointer:: node
  !.............................................................................
  if (self%temporary.or.self%allocated) return
  self%n = in%n
  self%l = in%l
  self%u = in%u
  self%ds = in%ds
  self%name = 'tmp'
  !$omp atomic
  nalloc = nalloc+1
  !$omp atomic
  id = id+1
  self%id = id
  if (present(label)) then
    if (verbose>0) print*,'  alloc: ', self%name, self%id, nalloc, ' '//label
  else
    if (verbose>0) print*,'  alloc: ', self%name, nalloc
  end if
  allocate (self%f(self%n(1), self%n(2), self%n(3)))
  self%temporary = .true.
  allocate (node)
  if (.not.associated(head)) head => node
  if (associated(tail)) tail%next => node
  node%scalar => self
  tail => node
END SUBROUTINE tempalloc

!> =============================================================================
SUBROUTINE tempdealloc (self, label)
  class(scalar_t):: self
  character(len=*), optional:: label
  !.............................................................................
  if (verbose>0) print*, self%id, self%temporary, self%allocated
  if (.not.self%temporary.or.self%allocated) return
  if (.not.associated(self%f)) return
  !$omp atomic
  nalloc = nalloc-1
  if (present(label)) then
    if (verbose>0) print*,'dealloc: ', self%name, self%id, nalloc, ' '//label
  else
    if (verbose>0) print*,'dealloc: ', self%name, self%id, nalloc
  end if
  deallocate (self%f)
  self%id = -self%id
  self%temporary = .false.
END SUBROUTINE tempdealloc

!> =============================================================================
SUBROUTINE scalar_cleanup
  type(node_t), pointer:: node, rm
  integer:: save
  !.............................................................................
  if (nalloc<=0) return
  !save = verbose
  !verbose = 1
  node => head
  do while (associated(node))
    call tempdealloc (node%scalar, 'cleanup')
    rm => node
    node => node%next
    deallocate (rm)
  end do
  nullify (head, tail)
  !verbose = save
END SUBROUTINE scalar_cleanup

!> =============================================================================
END MODULE scalar_mod
