# Phase 4: the complete object workspace, restored transactionally

Branch `agent/gui-phase4`, base `88c6ac21`. Implementation commit `c87b46da`; this
report follows in its own commit. Roadmap section 4.
## Changed behaviour

1. **Schema** (`TumultStudio2/session.py`, `SCHEMA_VERSION = 2`). `Workspace.from_dict`
   checks the version first, then every field type (no booleans as numbers, finite
   ascending ranges, equal-length finite knots), identities (unique spectrum/atom/result
   uids, positive unique mask-set uids, unique atom paths) and references (selections and
   the focus must name a saved object); a `SessionError` names the field. Spectra are
   selected by uid. `ObjectSnap`, the positional index, the `{"mset": …}` mask decoding
   and name-based relinking are gone; `MeasuredLines.from_dict` accepts only version 1.
2. **State persisted**: focus kind/uid, active tool, atom selection and tab, each wizard's
   step and jumped line, units, explicit main/nav x and y ranges, filter and display
   state. Restore routes the focus through `focus_spectrum`, `MaskTool.select_set`,
   `AtomTool.select` or `select_result`; a wizard, mask or atom panel is reactivated;
   range/filter/normalise gestures stay closed. A focus naming a skipped object is cleared.
3. **Prepare before replace** (`session.prepare`): spectra, mask sets, atoms and results
   are built with the workspace untouched; bad result records, mask sets or embedded data
   abort. Unreadable spectrum/atom files are skipped as `Skipped(kind, uid, name, path,
   reason)`; each wizard's `relink` clears references to them (reference observation, a
   balance's mask set; a correction's missing atom is noted and the result stays under
   Corrections). The footer lists counts, skipped items with path and reason and the
   number of cleared references; the log has one line per item.
4. **Commit** (`app.restore_workspace`): the live workspace is captured into a
   `PreparedWorkspace` (no file reads), the generation advances, non-detached jobs are
   cancelled, each wizard's `ResultJobs.current` is cleared, then `_apply_workspace`
   runs; if it raises, the captured workspace is re-applied and the footer says so.
   `ctx.restore_workspace(snapshot, name)` returns a `session.LoadReport`.
5. **Every collection replaced**, empty ones included, in the order spectra, masks,
   atoms, results, view/filter, focus, ranges; atom-only, mask-only, result-only, empty
   and all-spectra-missing slots are valid replacement states.
6. **Envelope shared**: `AbundanceWizard.snapshot/deserialise/relink/restore` with a
   per-tool `result_from_dict`; relinking stays explicit per tool. The Sessions list
   describes each slot (`session.summary`), a pre-phase one as "unsupported format".

## Schema summary
Top level: `schema_version`; `spectra` (uid, kind file/run/synth, name, path or embedded
data, meta, style, visible, conv, norm_knots, already_normalised); `spectrum_selected`;
`masks` (`sets` as `LineMaskSet.to_dict()` with baseline, `current`, `active`); `atoms`
(uid, path); `atom_selected`; `results` (`abundance`/`abundcorr`/`ionbal`, each
`{results, selected}`); `view` (xunit, yunit, main_x, main_y, nav_x, nav_y); `filter`;
`focus` (kind spectrum/linemask/atom/result/corr/ionbal, uid); `tools` (active, atom_tab,
steps, jumped). References are by identity only: `ref_uid` → spectrum uid,
`maskset_uid` → mask-set uid as a string, `atom_path` → atom path.

## Files
- `session.py` (schema, `prepare`, `LoadReport`, `summary`); `app.py` (`snapshot`,
  `prepare_workspace`, `_capture_live`, `_apply_workspace`, `restore_workspace`, slot
  list); `abundwizard.py`, `abundancetool.py`, `abundcorrtool.py`, `ionbalance.py`
  (item 6); `maskstore.py` (`restore(sets, current, active)`); `atomtool.py` (unused
  `to_dict` removed); `context.py`; `tumult/measured.py` (version check).
- `tests/test_session.py` (new, `.gitignore` negation): round trip, version rejection,
  18 malformed variants, summary, record version, preparation, missing files, broken
  embedded objects, empty/single-collection slots. `tests/test_ui_execution.py`: mixed
  round trip into a fresh document (text spectrum, `input/atoms/atom.ba06`), replacing an
  unrelated workspace and the empty slot, missing files and all-spectra-missing, four
  refused slots leaving state and generation intact, failed application rolled back, a
  worker finishing after a load. The legacy duplicate-name test went with the behaviour.
- `docs/GUIWorkflow.md`: "What a save slot holds", the load policy and the scope
  boundaries (forms stay global preferences; gestures reset; path-backed slots, no
  archives, hashes, export/import/rename/overwrite or resumable calculations).

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:
- `-m pytest tests/test_execution.py tests/test_ui_execution.py tests/test_store.py -q`:
  76 passed before; with `tests/test_session.py`, 94 passed after (13 session + 6 UI
  new, 1 legacy UI test removed), 33 s.
- `-m pytest tests/test_nlte_corrections.py tests/test_synthesis_contract.py -q`: 25 passed.
  `tests/test_linemask.py`: 9 failed / 3 passed / 2 skipped identically before and after
  (untracked `input/linemasks` data is absent). Serialization checks only; no physical
  validation, golden suites, build or browser.

## Limits and notes for later phases
- The main y range is a `DataRange1d`: the saved value is what the browser last synced
  and holds only until the next client autoscale; the nav y range is refitted on every
  nav refresh. Both are restored last, best effort. Slots saved before this phase are
  refused (version 1). `_capture_live` serialises the live results once, so a failed
  load of a large workspace costs more than a successful one.
- Phase 5: persist workflow references beside the workspace, not inside `Workspace`, with
  the `store` section/envelope fields and `CONFIG_VERSION` bump from the phase 3 notes.
  `restore_workspace` cancels non-detached jobs and advances the generation at commit, so
  a reattached handle needs a new `Job` under `ctx.generation` after a load;
  `LoadReport.cancelled_jobs` lists what the boundary cancelled.
- Phase 6: a request stored in `SpectrumSnap.meta` must stay JSON-serialisable; it is
  validated only as a dict. Phase 8: `AbundanceWizard.has_results` remains.
