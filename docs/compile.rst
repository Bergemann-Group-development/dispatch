Compiling
=========

To comile, go to the expriments/example directory, and use ``make infor``
to see macro names and values that would be used in ``make``::

   cd experiments/example
   make info

If the automatically selected compiler and options are OK, just do::

   make

If the compiler chosen by default is not appropriate, do instead, 
for example::

   make COMPILER=gfortran

.. toctree::
   :maxdepth: 3

   option_groups
   options

