AMR support
===========

The AMR actions on a task -- refine and derefine -- should only be
done by the thread that is about to update the task, since any other
alternative would lead to serious synchronization problems.  The
new level ``L-1`` tasks that may be needed to support a new level ``L`` task
can thus not be created until the next time the level ``L-2`` task that needs
to be refined is updated.

The ``refine_t`` data type has two methods related to checking for AMR
support (i.e., checking that all patches at level ``L`` have level ``L-1`` nbors
to get guard zone values from):

1. ``check_support()`` checks if a new (child) patch actually has support.
   If not, it issues a WARNING in the rank log file.  The task then takes
   guard zones values from ``L-2`` nbors, until the relevant ``L-1`` nbors have been
   created.
  
2. ``need_to_support()`` checks if a given ``L-2`` level patch needs to be
   refined, to provide support for existing level L patches.  If so, it
   creates such child patches(), along with any other refinement that may
   be detected a need for in ``refine_t%needed``.
  
Both methods use 2-byte integer maps to detect lack of support -- this
method will work unmodified for moving patches.
  
A status bit (``bits%support``) is available, and is used to communicate to
nbor patches about the need for support.

A new child patch is automatically sent vbor ranks when it is created, and
virtual patches are created there, as when any task arrives, with and id
that doesn't exist yet on the rank.  The ``bits%support`` informs the rank
that this is a new AMR child patch.

Any new task arriving to a rank generates a call to ``list_t%init_nbor_nbors()``,
which uses position information to first create an nbor list for the 
task itself, and then generates new nbor lists also for the nbors of the
task, since both sides of the nbor-nbor relation need to have consistent
nbor lists.

.. toctree::
   :maxdepth: 3

 
