
Initialization calls
====================
::

    experiment_t%init
      solver_t%init
        mhd_t%init
          self%idx%init        (obsolete!)
          self%initial%init    (obsolete!)
          timestep%init        [only stagger2]
          gpatch_t%init
            self%idx%init      (new)
            self%initial%init
            self%idx%init      (new)
            patch_t%init
              task_t%init
          force_t%init (obsolete!)
        extras_t%inig
          forces_t%init
            force_t%init
        validate%init

As much as possible should be inside framework files, avoiding
requiring that all ``$(SOLVER)%init`` contain a chain of specific calls.

We should thus consider moving calls to ``self%idx%init`` to ``gpatch_t``,
and doing it both before and after the call to ``self%initial%init``,
as is done in stagger2 (to pick up changes of the ``%mhd`` switch).

Any calls to ``self%initial%init`` and ``force%init`` in ``experiment_mod``
files, or in ``mhd_mod`` files, should be considered obsolete.

.. toctree::
   :maxdepth: 3

 
