# Phase 3: config transactions and save acknowledgements

Branch `agent/gui-phase3`, base `b052a1cc`. Implementation commit `1584029c`; this
report follows in its own commit. Roadmap section 3.

## Changed behaviour

1. **Missing versus broken config** (`TumultStudio2/store.py`). `_load` returns empty
   defaults only for `FileNotFoundError`. Any other read error, invalid JSON, a non-object
   document, a `version` other than `store.CONFIG_VERSION` (1), a section of the wrong
   type, a non-string recent file or a malformed/duplicate slot envelope raises
   `store.ConfigError`, whose message ends with the path and whose `.path` attribute
   holds it. Writers load inside their transaction, so a rejected read never becomes an
   empty replacement. New files are written with `"version": 1`; a file without the key
   is rejected as unsupported, with the message saying to move it aside. The slot
   `snapshot` content is not validated here (phase 4's workspace schema).
2. **Truthful acknowledgements** (`store.py`, `app.py` Sessions bridge). Write failures
   raise `ConfigError` with the OS cause (`Could not write config: No space left on
   device (<path>)`); the writer's temporary file is removed. `delete_session` returns
   whether the slot existed and does not rewrite an unchanged file. In `app.py`,
   `on_session_action` announces *Saved* only after `save_session` returned; a failed
   save, delete, load or refresh logs the traceback and shows `Save "name" failed: <cause>
   (<path>)` (or the delete/read equivalent) in the footer at danger level, leaving the
   previous file and the live canvas untouched for a retry. A delete of a slot that had
   already gone reports it at warn level.
3. **Serialised transactions.** `store.transaction()` (a context manager yielding the
   validated dict) takes a process-wide `threading.Lock` and an exclusive `flock` on
   `<config>.lock` (non-blocking retries for `LOCK_TIMEOUT_S` = 10 s, then `ConfigError`),
   loads, yields, and on normal exit replaces the file atomically from
   `<config>.<pid>.<thread>.tmp` if the dict changed. Every writer in `TumultStudio2/`
   already goes through the public store API (`set_pref`, `update_pref`, `set_section`,
   `add_recent`, `save_session`, `delete_session`: abundance/correction/ion-balance tools,
   synth setup/modal, grid manager, atom picker, fit tab, recent files), so all of them
   are now transactions; slot ids are allocated inside the lock. Readers take no lock.

## Files

- `TumultStudio2/store.py`: contract, `ConfigError`, `transaction`, writer-specific
  temporary files, `delete_session` return value.
- `TumultStudio2/app.py`: `render_sessions` and `on_session_action` only.
- `tests/test_store.py` (new; `.gitignore` negation added): missing versus eight broken
  configs, unreadable file, injected `os.replace` failure with temp-file cleanup, missing
  parent directory, delete outcome, unserialisable snapshot, 16 barrier-released threads,
  four interpreters writing 10 slots and preferences each against one path, lock timeout.
- `tests/test_ui_execution.py`: autouse fixture pointing every built document at a
  temporary config (the suite previously read and wrote the working directory's
  `.tumultstudio.json`); save acknowledgement with a failed then successful retry, delete
  failure, and an unsupported-version config blocking save/load/delete/refresh with the
  named cause, all through the real `session-actions` channel and `footer-status` bridge.
- `docs/GUIWorkflow.md`: new "Saved sessions and preferences" section.

## Verification

Worktree, `/Users/hoppe/Documents/dispatch/experiments/OpenMULTI3D/.venv/bin/python`
(`tumult` resolves into the worktree), `MPLCONFIGDIR=/tmp/tumult-mpl`:

- `-m pytest tests/test_execution.py tests/test_ui_execution.py -q`: 56 passed before;
  with `tests/test_store.py` added, 76 passed after (17 store + 3 UI new), 22 s.
- `pyflakes` on `store.py` and `test_store.py`: clean.
- Temporary directories, injected `os.replace`/permission failures and a Python
  subprocess writer only. No physical validation, golden suites, Fortran build or browser.

## Limits and notes for later phases

- The strict version check applies to existing user configs: a `.tumultstudio.json`
  written before this phase (the main checkout has one, with preferences and no slots)
  is rejected at `app.build()` until `"version": 1` is added or the file is moved aside;
  the error names the path. That is the roadmap's no-migration policy, so the failure is
  loud rather than silent.
- Store readers and writers outside the Sessions bridge (`store.add_recent` before a
  file/atom/run load, `set_pref`/`update_pref` from preference forms) now raise
  `ConfigError` instead of logging; those call sites were outside this phase's ownership,
  so such an error propagates to Bokeh's callback logging and aborts the action (the file
  is not opened). Deciding whether a recent-file failure should block opening a file
  belongs to whoever next touches those sites.
- `transaction()` is not reentrant; a store call inside another transaction deadlocks on
  the thread lock. Nested use does not exist in the tree.
- Phase 4: keep validating the snapshot in `session.py` and do it before any live change;
  `store.get_session` returns the raw dict or `None`, and `save_session(name, snapshot)`
  expects JSON-serialisable content (a `ConfigError` otherwise, written nothing). Use
  `store.transaction()` for any new multi-field update rather than several writers.
- Phase 5: workflow references belong in a new top-level section or inside the slot
  envelope; add the section to `_SECTIONS` (with its type) and the envelope fields to
  `_SLOT_FIELDS` so structural validation covers them, and bump `CONFIG_VERSION` when the
  layout changes. The lock covers local processes on one host; multi-host storage is out
  of scope.
