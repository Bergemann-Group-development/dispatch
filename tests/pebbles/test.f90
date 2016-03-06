PROGRAM test
  USE io_mod
  USE list_mod
  USE pebbles_mod
  implicit none
  class(link_t), pointer:: link
  class(pebble_patch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call io%init
  allocate (link, patch)
  link%task => patch
  call patch%test (link)
END PROGRAM test
