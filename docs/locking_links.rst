Neighbor list synchronization
-----------------------------

The link_t%nbor list is now changed nearly atomically, in that the new list is
first constructed separately, and is then connected with a single instruction,
pointing ``link%nbor`` to it.  The procedure ``list_t%init_nbors()`` creates 

1. a cached copy of the nbor list, which is used in ``list_t%send_to_vnbors``
2. a cached copy sorted by level, which is used in ``download_t%download_link()``

The nbor lists are also used in

3. ``list_t%check_ready()``, which checks those nbors that have ``nbor%needed``
   set, to see of their tasks are sufficiently advanced in time allow the linked
   task to be added to the ready queue
4. ``list_t%check_nbors()``, which checks those nbors that have ``nbor%is_needed``
   set, using ``check_ready()`` on each of them.  This is done after every active
   and virtual task update, since those are the only points in time when new
   tasks could become ready to update

The task link, which contains ``link%nbor``, must be protected by its OMP lock,
``link%lock``,
whenever the nbor list is updated or referenced.  The cached nbor lists are used
to keep these lock periods as short as possible
 

.. toctree::
   :maxdepth: 3

