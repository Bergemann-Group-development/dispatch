"""NLTE departures: solved, per grid node, or requested on demand.

Three objects state where the departure coefficients of a synthesis come from:

* :class:`Departures`: one equilibrium solve on one atmosphere as it stands on disk,
  together with what it was solved for (atmosphere, star, composition, model atom and
  wavelength quadrature). A synthesis given one reads the populations back and runs only
  the formal solution, for any window.
* :class:`DepartureGrid`: one :class:`Departures` per atmosphere node and microturbulence
  lattice point that a grid request interpolates between.
* :class:`DepartureRequest`: a description of equilibrium departures to be solved for
  whatever effective request consumes it. It holds the model atom and the quadrature
  policy but no request, so a star or abundance changed later changes the departures
  that are solved. Every abundance rung of a curve of growth then gets its own solve.

:func:`solve` solves eagerly and :func:`request` builds the deferred form::

    request = SynthRequest(config=config, atmosphere=atmosphere, star=star)
    solved = departures.solve(request, atom="input/atoms/atom.mg86g")
    blue = synthesis.synthesize(request.update(departures=solved, windows=(5160, 5190)))

    deferred = departures.request(request, atom="input/atoms/atom.mg86g")
    spectrum = synthesis.synthesize(request.update(departures=deferred, windows=(5160, 5190)))

The atom's bound-bound wavelength quadrature (:mod:`tumult.line_quadrature`) is prepared
once per operation, at the request's star and composition, on a MARCS atmosphere: the
atmosphere itself when it is one, else the MARCS node nearest the star in the grid or in
the atmosphere's ``preparation`` grid.
"""
from __future__ import annotations

import copy
import hashlib
from dataclasses import dataclass, field, replace
import json
from pathlib import Path
import time

MANIFEST = "departures.json"


def _atomic_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    tmp.replace(path)
    return path


@dataclass(frozen=True)
class Departures:
    """A solved departure run on disk, and what it was solved for.

    ``folder`` is the resolved run directory the spectrum phase reads back. ``name`` is the
    same run as the launcher names it, relative to ``datadir``, which is what a job must be
    given; only ``folder`` survives a change of working directory. ``atmosphere`` is the
    model the populations belong to, ``star`` the stellar parameters of the solve,
    ``composition`` its composition settings, ``atom`` the model atom it was requested for
    and ``quadrature`` the prepared wavelength quadrature (empty when the atom's own was
    used). ``settings`` is the full namelist of the solve.

    ``status`` says how this object was obtained in this process: ``"solved"`` by a native
    run, ``"reused"`` from the run cache or ``"loaded"`` from a saved record. It is not
    saved.
    """
    folder: Path
    name: Path
    atmosphere: str
    resolution: tuple = ()
    species: tuple = ()
    settings: dict = field(default_factory=dict, repr=False)
    composition: dict = field(default_factory=dict)
    convergence: dict = field(default_factory=dict)
    star: dict = field(default_factory=dict)
    atom: str = ""
    quadrature: dict = field(default_factory=dict, repr=False)
    preparation_atmosphere: str = ""
    status: str = field(default="loaded", compare=False)

    @property
    def saved(self):
        return (self.folder / "save").is_dir()

    @property
    def atom_settings(self):
        """``(atom_files, line_grid_file)`` the solve ran with."""
        atom = (self.settings or {}).get("atom_params", {}) or {}
        files = atom.get("atom_file") or []
        files = [files] if isinstance(files, str) else list(files)
        return [str(f) for f in files], atom.get("line_grid_file")

    def to_dict(self):
        return {"version": 2, "folder": str(self.folder), "name": str(self.name),
                "atmosphere": self.atmosphere, "resolution": list(self.resolution),
                "species": list(self.species), "settings": copy.deepcopy(self.settings),
                "composition": copy.deepcopy(self.composition),
                "convergence": copy.deepcopy(self.convergence), "star": dict(self.star),
                "atom": self.atom, "quadrature": copy.deepcopy(self.quadrature),
                "preparation_atmosphere": self.preparation_atmosphere}

    def save(self, path=None):
        """Persist the metadata beside the native saved populations, atomically."""
        return _atomic_write(Path(path) if path is not None else self.folder / MANIFEST,
                             self.to_dict())

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        if value.pop("version", None) != 2:
            raise ValueError("Unsupported departures record; open the run folder instead.")
        value["folder"] = Path(value["folder"]).expanduser().resolve()
        value["name"] = Path(value.get("name") or value["folder"].name)
        value["resolution"] = tuple(value.get("resolution") or ())
        value["species"] = tuple(value.get("species") or ())
        return cls(**value)

    @classmethod
    def load(cls, path):
        """Open saved departures from their record or their run folder.

        A run directory without a record is described from its resolved ``params.nml`` and
        native population files; its star and atom are then unknown.
        """
        path = Path(path).expanduser().resolve()
        manifest = path if path.is_file() else path / MANIFEST
        if manifest.is_file() and manifest.suffix == ".json":
            return cls.from_dict(json.loads(manifest.read_text()))
        folder = path.parent if path.is_file() else path
        params = folder / "params.nml"
        if not params.is_file():
            raise FileNotFoundError(f"No departures record or params.nml in {folder}")
        from .namelist import NameList
        settings = dict(NameList.from_file(params))
        atmos = settings.get("atmos_params", {}) or {}
        resolution = tuple(int(atmos[k]) for k in ("nx", "ny", "nz") if atmos.get(k))
        atom_files = (settings.get("atom_params", {}) or {}).get("atom_file") or []
        atom_files = [atom_files] if isinstance(atom_files, str) else list(atom_files)
        found = cls(folder=folder, name=Path(folder.name),
                    atmosphere=str(atmos.get("atmos_file") or ""), resolution=resolution,
                    species=atom_ids(folder), settings=settings,
                    composition=copy.deepcopy(settings.get("composition_params", {}) or {}),
                    atom=str(atom_files[0]) if atom_files else "")
        if not found.saved:
            raise ValueError(f"{folder} contains no native saved departures")
        return found


def atom_ids(folder):
    """Which atoms a run wrote departures for, from the files it wrote them in."""
    return tuple(sorted(p.name[len("pop_"):-len("_meta.txt")]
                        for p in (Path(folder) / "save").glob("pop_*_meta.txt")))


def _node_key_to_json(key):
    spatial, second = key
    return [list(spatial), second]


def _node_key_from_json(value):
    spatial, second = value
    return (tuple(float(v) for v in spatial), second)


@dataclass(frozen=True)
class DepartureGrid:
    """The departures of every node a grid request interpolates between.

    ``nodes`` pairs each node key of :func:`tumult.synthesis.select_nodes`
    (``((teff, logg, feh), vmic)``) with its :class:`Departures`; ``star`` is the stellar
    point the nodes were selected for and ``grid`` the atmosphere grid record.
    """
    star: dict
    grid: dict
    nodes: tuple
    atom: str = ""
    quadrature: dict = field(default_factory=dict, repr=False)
    preparation_atmosphere: str = ""

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(departures for _key, departures in self.nodes)

    @property
    def status(self):
        states = {departures.status for departures in self}
        return states.pop() if len(states) == 1 else "mixed"

    @property
    def atom_settings(self):
        return self.nodes[0][1].atom_settings

    def folders(self):
        """``{node key: departure folder}`` for the spectrum phase."""
        return {key: str(departures.folder) for key, departures in self.nodes}

    def to_dict(self):
        return {"version": 2, "star": dict(self.star), "grid": copy.deepcopy(self.grid),
                "nodes": [{"key": _node_key_to_json(key), "departures": dep.to_dict()}
                          for key, dep in self.nodes],
                "atom": self.atom, "quadrature": copy.deepcopy(self.quadrature),
                "preparation_atmosphere": self.preparation_atmosphere}

    def save(self, path):
        return _atomic_write(path, self.to_dict())

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        if value.pop("version", None) != 2:
            raise ValueError("Unsupported departure grid record.")
        nodes = tuple((_node_key_from_json(item["key"]),
                       Departures.from_dict(item["departures"]))
                      for item in value.pop("nodes"))
        return cls(nodes=nodes, **value)

    @classmethod
    def load(cls, path):
        return cls.from_dict(json.loads(Path(path).expanduser().read_text()))


@dataclass(frozen=True)
class DepartureRequest:
    """Equilibrium departures of ``atom``, solved for each request that consumes them.

    ``line_grid`` prepares the atom's bound-bound wavelength quadrature from its opacity
    profiles on the request's atmosphere (:func:`tumult.line_quadrature.prepare_line_quadrature`);
    off, the atom's own ``nq``/``qmax``/``q0`` quadrature is used. ``quadrature`` holds
    further generator settings (``R``, ``tol_dex``, ...).
    """
    atom: Path
    line_grid: bool = True
    quadrature: dict = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "atom", Path(self.atom).expanduser())
        object.__setattr__(self, "quadrature", dict(self.quadrature or {}))

    def prepare(self, request):
        """The prepared quadrature for ``request`` and the atmosphere it was prepared on.

        Returns ``(None, None)`` without ``line_grid``.
        """
        if not self.line_grid:
            return None, None
        from .line_quadrature import prepare_line_quadrature
        from .synthesis.operations import quadrature_atmosphere
        from .synthesis.workflows import preparation_root
        config, star = request.config, request.star
        atmosphere = quadrature_atmosphere(request)
        composition = request.abundances
        prepared = prepare_line_quadrature(
            preparation_root(config.to_prefs()), config.resolve(self.atom), atmosphere,
            vmic=star.vmic if star.vmic else 1.0, feh=star.feh,
            abundances=composition.override() if composition else "",
            relative_abundances=composition.relative if composition else True,
            **self.quadrature)
        return prepared, atmosphere

    def run_values(self, request, prepared=None):
        """The native run values selecting this atom, on ``prepared`` when given."""
        if prepared is None:
            return {"atoms": [str(request.config.resolve(self.atom))]}
        return {"atoms": [prepared.atom_file], "line_grid": prepared.grid_file}

    def to_dict(self):
        return {"atom": str(self.atom), "line_grid": self.line_grid,
                "quadrature": dict(self.quadrature)}

    @classmethod
    def from_dict(cls, value):
        return cls(**{k: v for k, v in dict(value).items() if k != "kind"})


KINDS = {"departures": Departures, "departure_grid": DepartureGrid,
         "departure_request": DepartureRequest}


def to_record(value):
    """A departures object as a tagged record, or ``None``."""
    if value is None:
        return None
    kind = {Departures: "departures", DepartureGrid: "departure_grid",
            DepartureRequest: "departure_request"}[type(value)]
    return {"kind": kind, **value.to_dict()}


def from_record(value):
    """The inverse of :func:`to_record`."""
    if value is None:
        return None
    value = dict(value)
    return KINDS[value.pop("kind")].from_dict(value)


def load(path):
    """Open saved :class:`Departures` (a record or run folder) or a :class:`DepartureGrid`."""
    path = Path(path).expanduser()
    if path.is_file() and path.suffix == ".json":
        value = json.loads(path.read_text())
        if "nodes" in value:
            return DepartureGrid.from_dict(value)
        return Departures.from_dict(value)
    return Departures.load(path)


def request(request, *, atom, line_grid=True, **quadrature):
    """Describe equilibrium departures of ``atom`` for consumers to solve on demand.

    ``request`` only validates the atom against the configuration; the returned
    :class:`DepartureRequest` does not keep it, so attaching it to an updated request
    solves at that request's star and abundances.
    """
    return _describe(request.config, atom, line_grid, quadrature)


def _describe(config, atom, line_grid, quadrature):
    path = config.resolve(atom)
    if not path.is_file():
        raise FileNotFoundError(f"Model atom not found: {path}")
    return DepartureRequest(atom, line_grid, quadrature)


def _status(folder, started):
    stamp = Path(folder) / "abund.txt"
    return "solved" if stamp.is_file() and stamp.stat().st_mtime >= started else "reused"


def solve(request, *, atom, line_grid=True, cancel_token=None, **quadrature):
    """Prepare ``atom``'s quadrature and solve its equilibrium departures for ``request``.

    An individual atmosphere or snapshot gives :class:`Departures`; a grid gives a
    :class:`DepartureGrid` with one solve per node and vmic lattice point the request's
    star interpolates between. Windows and the line list do not enter a departure solve.
    Both are saved beside their runs (``departures.json``; a grid's record lands in the
    cache under ``_departure_grids``).
    """
    from .synthesis.inputs import AtmosphereGrid
    from .synthesis.operations import native_request

    if request.departures is not None:
        raise ValueError("The request already carries departures; clear them to solve anew.")
    request.check()
    spec = _describe(request.config, atom, line_grid, quadrature)
    prepared, prep_atmosphere = spec.prepare(request)
    run_values = {"mode": "NLTE", "only_atom": True, **spec.run_values(request, prepared)}
    native = native_request(request.update(windows=()), **run_values)
    record = dict(atom=str(request.config.resolve(spec.atom)),
                  quadrature=prepared.to_dict() if prepared else {},
                  preparation_atmosphere=str(prep_atmosphere or ""))
    started = time.time() - 1.0
    if isinstance(request.atmosphere, AtmosphereGrid):
        return _solve_grid(request, native, record, started, cancel_token)

    from .synthesis.workflows import solve_request_departures
    solved = solve_request_departures(native, cancel_token=cancel_token)
    solved = replace(solved, star=dict(request.star), status=_status(solved.folder, started),
                     **record)
    solved.save()
    return solved


def _solve_grid(request, native, record, started, cancel_token):
    from .auxfuncs import update_dict
    from .synthesis.api import _dep_bucket, _dep_jobs, _save_provenance
    from .synthesis.atmos_select import select_nodes
    from .synthesis.executor import _run_batch
    from .synthesis.workflows import _prepare, resolve_request

    resolved = resolve_request(native)
    common, settings_list, _R = _prepare(native.target, native.prefs, [native.run_vals],
                                         native.grid, None, None, set_spectrum=False)
    settings = update_dict(settings_list[0], copy.deepcopy(native.expert_settings))
    nodes, _vmic = select_nodes(native.target, common["atmos_files"],
                                method=common["method"], vmic_step=common["vmic_step"],
                                absmet_files=common["absmet_files"])
    jobs = _dep_jobs(settings, nodes, _dep_bucket(Path(common["top"]), settings))
    runs = _run_batch(jobs, common["launch_kwargs"], common["max_threads"],
                      cancel_token=cancel_token)
    failed = [key for key, _path, _nml in jobs if runs.get(key) is None]
    if failed:
        raise RuntimeError(f"Departure runs failed for atmosphere nodes: {failed}")
    _save_provenance(runs, nodes, settings)
    root = Path(resolved.execution_root)
    solved = []
    for (key, teff, logg, feh, vmic, _weight, atmos_path, _absmet), (_k, name, nml) in \
            zip(nodes, jobs):
        folder = (root / runs[key].folder).resolve()
        node = Departures(
            folder=folder, name=Path(name), atmosphere=str(atmos_path), species=atom_ids(folder),
            settings=nml, composition=copy.deepcopy(nml.get("composition_params", {}) or {}),
            star={"teff": teff, "logg": logg, "feh": feh, "vmic": vmic},
            status=_status(folder, started), **record)
        node.save()
        solved.append((key, node))
    grid = DepartureGrid(star=dict(request.star), grid=native.grid, nodes=tuple(solved),
                         **record)
    digest = hashlib.sha1(json.dumps([grid.star, sorted(grid.folders().values())],
                                    sort_keys=True).encode()).hexdigest()[:12]
    grid.save(Path(resolved.output_root) / "_departure_grids" / f"{digest}.json")
    return grid
