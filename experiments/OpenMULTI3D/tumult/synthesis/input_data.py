"""The published input data that does not live in git: the Keeper share folders.

Besides the atmosphere grids (:mod:`.grid_catalog`), OpenMULTI3D's larger inputs are
folders of plain files on one public Seafile share (MPCDF Keeper): the master line lists,
the Gaia-ESO NLTE line list, the public model atoms, the metal-opacity (absmet) tables,
the line masks and observed spectra (the KPNO solar flux atlas, the IAG disc-resolved solar
spectra and the Gaia FGK benchmark stars v3 library) and window files with the wavelength
coverage of spectroscopic surveys. Each :class:`Dataset` names a share folder and where its files belong
under the checkout's ``input/`` directory; the line masks keep their per-element
sub-folders. The listing, planning and resumable transfer are those of
:mod:`.linelist_catalog`, which also serve the TumultStudio Preferences pages.

The command line and the Python entry point are :mod:`tumult.downloads`, which serves
these datasets and the atmosphere grids of :mod:`.grid_catalog` alike.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..launch import CHECKOUT
from . import linelist_catalog
from .linelist_catalog import ShareFolder

_KEEPER = dict(server=linelist_catalog.SOURCE.server, token=linelist_catalog.SOURCE.token,
               citation=linelist_catalog.SOURCE.citation)


@dataclass(frozen=True)
class Dataset:
    """A share folder and its default place under ``input/``."""
    id: str
    source: ShareFolder
    folder: str            # relative to the checkout's input/ directory

    @property
    def default_dest(self):
        return CHECKOUT / "input" / self.folder


DATASETS = {d.id: d for d in (
    Dataset("linelists", linelist_catalog.SOURCE, "linelists/master_linelists"),
    Dataset("ges", ShareFolder(id="keeper-ges", name="Gaia-ESO NLTE line list",
                               path="/linelists/GES/", **_KEEPER), "linelists"),
    Dataset("atoms", ShareFolder(id="keeper-atoms", name="Model atoms",
                                 path="/atoms/public/", **_KEEPER), "atoms/public"),
    Dataset("absmet", ShareFolder(id="keeper-absmet", name="Metal opacity tables",
                                  path="/absmet_files/", **_KEEPER), "absmet_files"),
    Dataset("linemasks", ShareFolder(id="keeper-linemasks", name="Line masks",
                                     path="/linemasks/", nested=True, **_KEEPER),
            "linemasks"),
    Dataset("windows", ShareFolder(id="keeper-windows", name="Survey window files",
                                   path="/windows/", **_KEEPER), "windows"),
    Dataset("kpno", ShareFolder(id="keeper-kpno", name="KPNO solar flux atlas",
                                path="/spectra/KPNO/", **_KEEPER), "spectra"),
    Dataset("iag", ShareFolder(id="keeper-iag", name="IAG disc-resolved solar spectra",
                               path="/spectra/IAG/", **_KEEPER), "spectra/IAG"),
    Dataset("gbs", ShareFolder(id="keeper-gbs", name="Gaia FGK benchmark stars v3",
                               path="/spectra/GBSv3/", **_KEEPER), "spectra/GBSv3"),
)}
