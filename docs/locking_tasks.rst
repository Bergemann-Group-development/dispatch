Task mem array synchronization
------------------------------

During a task update, the only memory slot that is written to is task_t%new,
while potentially all memory slots are needed to produce the new values.

Before a task update, all memory slots (1:nt-1) of nbor (source) tasks may in 
principle be needed, when computing guard zone values for the (target) task that
the thread is going to update.

If the nbor task is dormant (with bits%busy not set), it could possibly start
computing values after its time slot info has been acquired, but since the 
download for a single task takes less than 1% of the update, there is no chance
whatsoever that the task has time to both finish one update, and start producing
new mem array values for the next time slot (number 1 in terms of the iit array).

If the task is active (with bits%busy set), and is just about to finish its
update, it might have time to do that, changing in the process its time slot
indices and its task_t%time, but with no chance to finish one more update.

Hence, after retrieving (atomically) the time slot information from a source task,
using the task_t%timeslots() procedure, the thread working on computing guard
zone values for the target task can be sure that it can access slots (1:nt-1)
without conflicting with updates to the source task.  

The only possible exception would be if the task was suspended in the middle
of the download_t%different() or download_t%same() procedure call, and did not
wake up until the source task had updated several times.   This might possibly
(but still with very low probability) happen if one uses hyperthreading with a
large factor (such as using 20 threads on 2 cores).

It is important, however, to make sure that the source%t and souce%iit values
are retrieved only once, and that they are consistent, but if the local iit()
values are based on a source%it that is updated as the last thing in task_t%rotate(),
then even if the slot information is being changed by some thread execting 
source%rotate() during the target%timeslots() call, the target will still
receive consistent (1:nt-1) values.  

Consider the situation before the source update finishes:::

  iit(1) iit(2) iit(3) iit(4) iit(5)
    4      5      1      2      3
                       (it)   (new)
   t(1)   t(2)   t(3)   t(4)   t(5)
   0.3    0.4    0.0    0.1    0.2
   (1)    (it)  (new)   (4)    (5)
 
All that happens when the source update finishes is that t(3) is set to 0.5,
the (new) mem slot is filled with updated variable values, and that the iit
array is shifted left, so it reads::
 
  iit(1) iit(2) iit(3) iit(4) iit(5)
    5      1      2      3      4
                       (it)   (new)

The shift of iit(:) may be replaced by knowing only the value of "it", which
is set atomically.   So, if "it" is incremented (cyclically) right after the
time slot that it corresponds to has been updated, then a thread downloading
to a target gets a set of times to interpolate between that is consistent with
the content of the corresponding memory slots.  Should the thread pick up "it"
a microsecond too early it still gets useful data, since already the check_ready
test that happened typically a long time ago made the judgement that the task
was ready to be updated.

.. toctree::
   :maxdepth: 3

