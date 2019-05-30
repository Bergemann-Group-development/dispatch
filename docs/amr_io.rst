AMR save/restore
=======================================

A saved AMR snapshot must contain -- effectively -- a list of tasks to read in,
initialize, and fill with values from the snapshot.

In the case of for example the ``pdd_zoom`` experiment, the task list consists of
more than one piece, with different types of task in each piece.  Sufficient info 
about this must be carried by the snapshot meta-data.

In general, a "component" that is part of the simulation setup ("scene") may need
to have its own save/restore (output/input).  If so, it also means that patches
and tasks need to be identified as belonging to a component, either by being part
of a partial task list, or by having an identifier string (or integer).  The former
is not convenient, since on the one hand we want to have all tasks in a rank-global
task list, and on the other hand we need to be able to add and remove tasks in only
one place, not in the global and partial lists separately.

One could possibly rely on that tasks needing special I/O treatment probably also
need extra, non-generic parameters, and hence need to belong to a certain type of
task, which then would naturally have its own input/output procedure.

In any case, when reading in AMR tasks, one needs to know the task type to be
able to create the task, into which the data is to be read.  The sequence of steps
needs to be

1. read the meta-data of a task
2. create the task
3. load data into the task
    
Clearly, the meta-data needs to contain information both about the task type, and
about the location of the task data on disk.

It would be nice if this could be semi-automatic, so detection of task type was
done by a required data type method.

A part of this issue is also the arrangement of different task types.  We should
be able to allow different solvers to exist, without insisting that one must be
the extension of another.  They should all be extensions of task_t, and some should
all be extensions of patch_t:::


                 task_t
                 /  |  \
                /   |   \
         task1_t    |    task2_t
                    |
                 patch_t
                /   |   \
               /    |    \
       patch1_t     |     patch2_t

.. toctree::
   :maxdepth: 3

 
