PROGRAM amr_example
  USE io, only: init_io
  USE amr_mod, only: amr
  USE refine_mod
  USE load_balance_mod
  implicit none
  type(amr_t), pointer:: a
!...............................................................................
  do_debug = .true.             ! debug
  mpi_size = 5                  ! fake number of ranks
  call init_io                  ! for tracing and debug
  a%maxamr = 30
  call init_amr                 ! allocate 30 blocks (we need 25)
  a => amr

  call make_sons(a, 1)          ! create 8 sons of the root block => level 2
  call print_amr (a)            ! print status

  call change_rank(a, 5, 1)     ! change the rank of the 4th son to 1
  call print_amr (a)            ! print status

  a%patch(5)%refine = .true.    ! mark block 2 for refinemnt => level 3
  call refine (a)               ! carry out refinement
  call print_amr (a)            ! print status

  call make_sons(a, 11)         ! refine one of the level 3 blocks => level 4
  call print_amr (a)            ! print status

  call load_balance (a)         ! load balance
  call print_amr (a)            ! print status

  call deallocate_amr (a)
END PROGRAM amr_example
