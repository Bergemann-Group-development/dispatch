Python ``dispatch`` module
==========================

Try for example::

  import dispatch
  help(dispatch)
  help(dispatch.snapshot)
  s=dispatch.snapshot(<TAB>

To read the 1st snapshot from ``data/``, ``data/run``, and ``../data/run``, do::

  s=dispatch.snapshot(1)
  s=dispatch.snapshot(1,run='run')
  s=dispatch.snapshot(1,'run',data='../data')

The resulting ``s`` is an instance of the snapshot class.  To see what it
offers, do::

  dir(s)

