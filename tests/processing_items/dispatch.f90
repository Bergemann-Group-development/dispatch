PROGRAM dispatch
  USE item_mod
  USE task_mod
  USE experiment_mod
  implicit none
  class(experiment_list_t), pointer:: task_list
  class(item_t), pointer:: item
  !.............................................................................
  allocate (task_list)
  call task_list%init ('my experiment')
  do while (task_list%size > 0)
    item => task_list%ready_queue%head%item
    select type (item)
    type is (item_t)
      call item%update
    type is (task_t)
      call item%update
    type is (experiment_t)
      call item%update
    end select
    call task_list%check (item)              ! check if item finished, etc
  end do
END PROGRAM
