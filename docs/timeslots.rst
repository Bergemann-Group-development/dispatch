Timeslots
==========

Each task has associated with it two sets of information that are also accessed
by other threads:

1) the ``task%t(:)``, ``task%dt(:)``, ``task%time``, and ``task%it``
   (equivalently ``task%iit``) info about time slots

2) the ``link%nbor`` nbor lists, which lists the nbor tasks that the task depends on

The 1st set of values are accessed twice:  First in ``list_t%check_nbors()`` and
``check_ready()``, which uses task_t%time and task_t%dtime to determine if the nbors
of a task are sufficiently advanced in time to consider the task ready to update.

Later, when the guard zone values are downloaded, all of the first ``(nt-1)`` values
of ``task_t%t(:)`` and ``task_t%dt(:)`` are needed.

The 2nd set of values (or pointers) are used in both of these steps as well.
However, since the nbor lists are not changed, or else are cached, there is
no problem or need for locking in this context, except very briefly, while
making the sorted cache copy and, correspondingly, when changing the nbor
lists in connection with refine / derefine.
