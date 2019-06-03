Garbage collection
==================

When tasks are being added or removed it is convenient if threads the are still
assuming that the tasks exist can finish their work, without errors or task
locking.  This can be achieved by maintaining a count ``task_t%n_needed``, which
keeps track of how many threads that currently have a lock at an nbor list where
the task is a member.

If the number is larger than zero, a task that should be deleted is instead added
to a list of garbage task, from which it is removed and deleted when the count
reaches zero.

There are only three procedure that create or delete nbor lists: ``copy_nbor_list``,
``sort_nbors_by_level``, and ``remove_nbor_list``.   These maintain the ``n_neeeded``
count, using atomic increments and decrements.
