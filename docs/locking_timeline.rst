Locking timeline
----------------

Schematic timeline and call hierarchy showing lock set & unset::

  task_t%update
    refine_t%check_current
      refine_needed
        need_to_support
          tlist%lock%set                                                         +|
          link%lock%set                                                 +|        |
          ..                                                             |        |
          link%lock%unset                                               -|        |
          tlist%lock%unset                                                       -|
      selective_refine
        make_child_patch
          parent%lock%set                                               +|
          parent%lock%unset                                             -|
          child%link%init
          child%init
            download_link
          tlist%lock%set                                                         +|
          tlist%append_link                                                       |
          tlist%init_nbors (parent%link)                                          |
            link%lock%set                                                         .
            do                                                                    .
              link%add_nbor_by_rank                                               .
              link%set_status_and_flags
            link%copy_nbor_list (... nbors)
            link%remove_nbor_list (... %nbors_by_level)
            link%sort_nbors_by_level (... %nbors_by_level)
            link%increment_needed (... %nbors_by_level)
            link%decrement_needed (nbors)
            link%remove_nbor_list (old_head)
          init_nbors (child_link)
             ditto
          init_nbor_nbors (child_link)
             link%lock%set
             self%init_nbors (link)                                     +|
             copy_nbor_list (link%nbor, nbors)                           |
             link%lock%unset                                            -|
             for nbor in link%nbor
               init_nbors (nbor%link)
             link%remove_nbor_list (nbors)
          reset_status
          check_support()
            tlist%lock%set
            do link in tlist
              ...
            tlist%lock%unset
          send_to_vnbors
          check_nbor_nbors
            for nbor 
              check_nbors (nbor%link)                                             .
                ...                                                               .
                check_ready                                                       .
            check_ready (link)                                                    |
          tlist%lock%unset                                                       -|
        remove_patch
      check_support(tlist,link)
    task%dnload
      download_link
        link%lock%set
        link%copy_nbor_list (nbors)
        link%increment_needed (nbors)
        link%lock%unset
        do
          source%lock%set
          same
          different
          source%lock%unset
        end do
        link%decrement_needed (nbors)
        link%remove_nbor_list (nbors)
    task%update
    task%rotate
    task%info
    send_to_vnbors
    check_nbors
      link%lock%set
      copy_nbor_list (... nbors)
      increment_needed (nbors)
      link%lock%unset
      for nbor ..
        check_ready (nbor%link)
          link%lock%set
          do .. nbor
            ..
          link%lock%unset
        queue_by_time (link)
          tlist%lock%set
    task%has_finished

.. toctree::
   :maxdepth: 3

