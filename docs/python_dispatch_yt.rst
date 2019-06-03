Python ``dispatch.yt`` module
=========================================

``dispatch.yt`` is an interface to the `YT Project <https://yt-project.org/>`_ Python package.

Try for example::

  import yt
  import dispatch.yt

To read the 1st snapshot from the ``../data/run``, do::

  help(dispatch.yt.snapshot)
  ds=dispatch.yt.snapshot(1,run='run',data='../data')

The resulting ``ds`` is an instance of the YT data set class.
For documentation of the various YT commands below, use ``help(command)``.
To make a slice plot, try::

  fld='density'
  slc=yt.SlicePlot(ds,fld=fld,axis=2,center=[.5.,.5,.5])
  slc.set_log(fld,True)
  slc.set_colorbar_label(fld,'Density')
  slc.set_xlabel('')
  slc.set_ylabel('')
  slc.annotate_grids(edgecolors='white',draw_ids=False)
  slc.show()
  slc.save('slice_plot.png')

