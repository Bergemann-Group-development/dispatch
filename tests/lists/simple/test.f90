!==============================================================================
PROGRAM test
 USE sll_mod
 USE dll_mod
 type(sll_t):: sll
 type(dll_t):: dll
 !------------------------------------------------------------------------------
 call sll%test
 call dll%test
END PROGRAM test
