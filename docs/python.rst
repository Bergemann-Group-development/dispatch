Python support
==============

To load Python modulesi, add the DISPATCH directory ``utilities/python/``
to the Python path (environment variable ``PYTHONPATH``), and do::

    import dispatch
    import dispatch.yt
    import dispatch.select
    import dispatch.graphics

The modules are self-documented, cf:

.. toctree::
   :maxdepth: 3

   python_dispatch.rst
   python_dispatch_yt.rst
   python_dispatch_select.rst
   python_dispatch_graphics.rst

The source code may be found in ``utilities/python/dispatch/``, and is
also available via the :ref:`auto-generated documentation <doxy>` (the
"Files" menu contains source code).
