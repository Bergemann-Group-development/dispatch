
Task locking
============

DISPATCH uses OpenMP nested locks, with the API defined in ``omp/omp_locks_mod.f90``.

A basic rule that needs to be respected with locks is to avoid that two objects 
of the same class are locked at the same time by a task -- the situation below
can clearly lead to a deadlock ::

  thread 1:
    lock A(1)
      lock & unlock A(2)
    unlock A(1)
  thread 2:
    lock A(2)
      lock & unlock A(1)
    unlock A(2)
    
On the other hand, if a thread locks one type of lock, and then a whole set
of another type of locks inside the first one, this cannot lead to a deadlock,
since either the outer locks are the same, and only one thread can do the set,
or they are different, and they can negotiate.  Hence locking the task list 
avoids any deadlock that could potentially be triggered if threads were not
locked out.

.. toctree::
   :maxdepth: 3

   locking_tasks
   locking_links
   locking_timeline
   locking_others
   locking_guard_zones
   locking_measures
   locking_nbor_lists
   locking_task_memory
   locking_timings
