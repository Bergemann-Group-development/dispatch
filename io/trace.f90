SUBROUTINE trace_begin (label)
  USE trace_mod
  implicit none
  character(*):: label
  call trace%begin (label)
END SUBROUTINE trace_begin

SUBROUTINE trace_end ()
  USE trace_mod
  implicit none
  call trace%end ()
END SUBROUTINE trace_end
