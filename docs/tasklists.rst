Task lists
===========

Tasks lists are fundamental to DISPATCH, and the list nodes (data type link_t)
contain several types of pointers that define relations between tasks; e.g. a
subset with time order, or a subset of active tasks.

The task list should not be required to know any specifics of tasks, which are
referred to with pointers of type experiment_t; the anonymous top tier task,
which may be anything between the top of a hierarchy with many layers (e.g.
task_t -> patch_t -> gpatch_t -> mhd_t -> rmhd_t -> solver_t), or just two
layers (task_t -> experiment_t).

The task_list_t data type extends the list_t data type with procedures that can
handle tasks in the context of task lists; e.g. specifying what happens in the
context of unpacking a task MPI package from another MPI process.

Procedures
----------

As is, the task_list_t data type contains a number of procedures that have to
do with message sending and receiving.  These could possibly be split off into
task_mesg_t data type.  Alternatively, the procedures that in effect implement
dispatcher method=0 could be split off into a dispatcher_method0_t data type,
leaving task_list_t to effectively be the task_mesg_t.

The list_t data type contains (or should contain only) procedures that manipulate
lists of link_t data types; inserting or deleting nodes in lists, etc.

As is, the list_t data type also contains procedures that rely on tasks having
meshes, e.g. for constructing neighbor lists.  These could perhaps with advantage
be split off into a patch_list_t data type.

tasks, one should be able to define task type specific relations that define 
when a source task is ahead of a target task.   This requires the functions that
call is_ahead_of are at a level where they are aware of experiment_t and all
sub-levels.::

|    dispatcher_t                                                  dispatcherX_mod
|      task_list_t check_ready                                     task_list_mod
|        experiment_t is_ahead_of refine                           experiment_mod
|          solver_t                                                solver_mod
|            gpatch_t                                              gpatch_mod
|              patch_t                                             patch_mod
|                task_t                                            task_mod

For task refinement, a similar desire exists.  The procedure that defines if a 
task should be "refined" (whatever that means) should be aware of all levels of
the hierarchy.  Or else, one should be able to overload the refine procedure 
itself, at any level.

Refinement
----------

How should refinement procedures and the dispatchers really work together?
Currently, the task list updater calls are refine procedure, as part of the 
task_list_t%update procedure, but this is awkward, since it means that one
allows a procedure that really deals with a single link in a task list to 
affect the task list itself.  It would be safer and more consistent if the 
level that loops over task list events (i.e., the dispatcher level) is the
one that also considers refinement.

But if the refine_mo should be able to manipulate the task list with dispather0,
it needs to know about task list, which creates a Catch 22, since it is also 
called from inside task_list_mod.::

|    dispatcher_t                                                  dispatcherX_mod
|    task_list_t                                                   task_list_mod
|    list_t                             check_ready                list_mod
|                   experiment_t                                   experiment_mod
|                   rt_solver_t                                    rt_solver_mod
|                   rt_t               (is_ahead_of)               rt_mod
|    refine_t ---------------------------------------------------- refine_mod
|                   solver_t              |                        solver_mod
|                   mhd_t                 |                        mhd_mod
|                   gpatch_t              |                        gpatch_mod
|                   patch_t               |                        patch_mod
|    link_t                               |                        link_mod
|                   task_t              is_ahead_of                task_mod

.. toctree::
   :maxdepth: 3

 
