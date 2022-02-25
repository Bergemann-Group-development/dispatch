!===============================================================================
!> Module used in gpatch_mod (only), to connect arrays at levels above
!===============================================================================
MODULE connect_mod
  type connect_t
    real, pointer, contiguous, dimension(:,:,:)     :: ex  => null()
    real, pointer, contiguous, dimension(:,:,:)     :: ey  => null()
    real, pointer, contiguous, dimension(:,:,:)     :: ez  => null()
    real, pointer, contiguous, dimension(:,:,:,:)   :: em  => null()
    real, pointer, contiguous, dimension(:,:,:,:)   :: uin => null()
    real, pointer, contiguous, dimension(:,:,:,:)   :: q   => null()
    real, pointer,             dimension(:,:,:,:)   :: bf  => null()
    real, pointer, contiguous, dimension(:,:,:,:,:) :: fl  => null()
  end type
END MODULE connect_mod
