// TumultStudio workspace: the shell's own behaviour — the Objects and Processes sidebar,
// the Open / New / Edit menus and their import dialogs, saved workspaces, the inspector's
// Edit tab, the spectrum tools (range, filter, normalise, line list, line masks) and the
// draggable bottom panels. Binds after bridge.js; the analysis and dialog scripts bind
// after this one, so a click that both the object list and a tool act on reaches the
// object list first.
(function () {
  var TS = window.TumultStudio;
  var model = TS.model, wireModal = TS.wireModal, seedResults = TS.seedResults;
  var wireNav = TS.wireNav, wireMenu = TS.wireMenu, ewAction = TS.ewAction, busySubmit = TS.busySubmit;

  wireModal(null, null, 'workspace-load-modal');
  wireModal(null, null, 'job-modal');
  wireModal(null, null, 'recovery-modal');

  // a Processes row names when its job started (data-since, epoch seconds); after 10 s
  // the elapsed time follows its state, refreshed here rather than by the server
  function formatElapsed(seconds) {
    seconds = Math.floor(seconds);
    if (seconds < 60) return seconds + ' s';
    if (seconds < 3600) return Math.floor(seconds / 60) + ' min';
    var minutes = Math.floor(seconds % 3600 / 60);
    return Math.floor(seconds / 3600) + ' h ' + (minutes < 10 ? '0' : '') + minutes + ' min';
  }
  setInterval(function () {
    var now = Date.now() / 1000;
    document.querySelectorAll('#processes .process-elapsed[data-since]').forEach(function (el) {
      var seconds = now - parseFloat(el.dataset.since);
      el.textContent = seconds >= 10 ? ' · ' + formatElapsed(seconds) : '';
    });
  }, 1000);

  // Edit tab: "Fit to target" / "Select full range" (rendered dynamically) ->
  // "edit-actions"; the Info tab's "Rerun from this…" -> "object-actions"
  document.addEventListener('click', function (ev) {
    var rerun = ev.target.closest('.object-rerun');
    if (rerun) {
      var o = model('object-actions');
      if (o) o.setv({ value: 'rerun:' + rerun.dataset.uid + ':' + Date.now() });
      return;
    }
    var kind = ev.target.closest('.edit-fit-btn') ? 'fit' :
               ev.target.closest('.edit-mask-btn') ? 'mask' : null;
    if (!kind) return;
    var m = model('edit-actions');
    if (m) m.setv({ value: kind + ' ' + Date.now() });
  });

  // the import modals; all are opened from the Objects "Open" menu. Each seeds its
  // picker's listing on open (the results container names its bridge)
  var openSpectrumModal = wireModal(null, 'close-file', 'file-modal', function () {
    seedResults('file-results');
  });
  var openLinemaskModal = wireModal(null, 'close-linemask', 'linemask-modal', function () {
    seedResults('linemask-results');
  });
  var openAtomModal = wireModal(null, 'close-atom', 'atom-modal', function () {
    seedResults('atom-file-results');
  });
  wireNav('#linemask-modal');

  function openMaskSource(source) {
    var modal = document.getElementById('linemask-modal');
    modal.classList.toggle('mask-open', source === 'file');
    modal.classList.toggle('mask-new', source === 'blank');
    modal.querySelector('.modal-title').textContent = source === 'file' ? 'Open line mask set' : 'New line mask set';
    modal.querySelector('[data-page="' + source + '"].settings-nav-item').click();
    openLinemaskModal();
  }
  wireMenu('open-file', 'open-file-menu', function(item) {
    if (item.dataset.import === 'linemask') openMaskSource('file');
    else if (item.dataset.import === 'atom') openAtomModal();
    else openSpectrumModal();
  });
  wireMenu('new-object', 'new-object-menu', function(item) {
    if (item.dataset.new === 'mask') openMaskSource('blank');
    else ewAction('synth-new');
  });
  // Edit: Rename… / Duplicate / Delete on the focused object; the items carry its kind, uid
  // and name, and are rendered again as the menu opens (editmenu.py)
  function editAction(op, item, name) {
    var m = model('edit-menu-actions');
    if (m) m.setv({ value: [op, item.dataset.kind || '', item.dataset.uid || '', name || '',
                            Date.now()].join('\u0000') });
  }
  wireMenu('edit-object', 'edit-object-menu', function(item) {
    var op = item.dataset.edit;
    if (op === 'rename') {
      var name = window.prompt('Rename “' + item.dataset.name + '” to:', item.dataset.name);
      if (name !== null) editAction(op, item, name);
    } else editAction(op, item);
  }, function() { editAction('refresh', {dataset: {}}); });
  document.getElementById('line-list-btn').addEventListener('click', function() { ewAction('lineopen'); });
  // session actions (save/load/delete/refresh) go through the "session-actions"
  // channel: space-delimited so a save name may contain spaces or colons
  function sessionAction(kind, arg) {
    var m = model('session-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }

  wireModal('open-session', 'close-session', 'session-modal', function () {
    sessionAction('refresh');  // re-read saves (another tab may have added some)
  });

  // load the currently-open directory as a tumult run
  document.getElementById('load-run').addEventListener('click', function () {
    var m = model('object-actions');
    if (m) m.setv({ value: 'loadrun::' + Date.now() });
  });

  // add a new named save of the current workspace
  function addSave() {
    var inp = document.getElementById('session-name');
    sessionAction('save', (inp.value || '').trim());
    inp.value = '';
  }
  document.getElementById('session-save').addEventListener('click', addSave);
  document.getElementById('session-name').addEventListener('keydown', function (e) {
    if (e.key === 'Enter') { e.preventDefault(); addSave(); }
  });

  // load / delete an existing save slot — delegated so it works in both the Sessions
  // modal and the inline empty-state loader shown over the canvas
  document.addEventListener('click', function (ev) {
    var del = ev.target.closest('.session-del');
    if (del) { ev.stopPropagation(); if (confirm('Delete this saved workspace?')) sessionAction('delete', del.dataset.sid); return; }
    var item = ev.target.closest('.session-item');
    if (item) {
      var panel = item.closest('.empty-card, dialog');
      panel.querySelectorAll('.session-item.selected').forEach(function(r) { r.classList.remove('selected'); });
      item.classList.add('selected'); panel.querySelector('.session-load').disabled = false;
    }
    var load = ev.target.closest('.session-load');
    if (load) {
      var selected = load.closest('.empty-card, dialog').querySelector('.session-item.selected');
      if (selected) sessionAction('load', selected.dataset.sid);
    }
    var choice = ev.target.closest('[data-workspace-load]');
    if (choice) sessionAction(choice.dataset.workspaceLoad);
  });

  // empty-state quick actions: reuse the existing Open spectrum modal / Synth toggle
  document.getElementById('empty-open').addEventListener('click', openSpectrumModal);
  document.getElementById('empty-synth').addEventListener('click', function () { ewAction('synth-new'); });
  document.getElementById('empty-atom').addEventListener('click', openAtomModal);

  // where a clicked line-list row sat when it was clicked ({lam0, off} in pixels from the
  // top of the scroll box), set on the click and consumed by the render it provoked
  var linelistAnchor = null;

  // delegated clicks on the object list: buttons hide/remove, the row selects
  document.addEventListener('click', function (ev) {
    var m = model('object-actions');
    if (!m) return;
    // the abundance tables share the .line-sort header class but report through their own
    // channels (analysis.js), so ignore their headers here
    var sortth = ev.target.closest('.line-sort');
    if (sortth && !sortth.closest('.abund-table') && !sortth.closest('.atom-table')
        && !sortth.closest('.grid-models')) {
      m.setv({ value: 'linesort:' + sortth.dataset.sortcol + ':' + Date.now() });
      return;
    }
    var row = ev.target.closest('.line-row');
    if (row) {
      // the jump re-cuts the table to a narrow window around this line, so remember where
      // the row sits now: the re-render puts it back there (see linelistAnchor below)
      var lslot = document.getElementById('linelist-table');
      linelistAnchor = lslot ? { lam0: row.dataset.lam0, ion: row.dataset.ion,
                                 off: row.getBoundingClientRect().top
                                      - lslot.getBoundingClientRect().top } : null;
      m.setv({ value: 'linejump:' + row.dataset.lam0 + '|'
                      + encodeURIComponent(row.dataset.ion) + ':' + Date.now() });
      return;
    }
    var swatch = ev.target.closest('.color-swatch');
    if (swatch) {
      m.setv({ value: 'color:' + swatch.dataset.color + ':' + Date.now() });
      return;
    }
    var undo = ev.target.closest('[data-undo]');
    if (undo) {
      // a removal's Undo in the status message (removal.py)
      var r = model('removal-actions');
      if (r) r.setv({ value: 'undo ' + undo.dataset.undo + ' ' + Date.now() });
      return;
    }
    var copy = ev.target.closest('.job-copy');
    if (copy) {
      if (navigator.clipboard) navigator.clipboard.writeText(copy.dataset.copy);
      return;
    }
    var jobAction = ev.target.closest('.job-action');
    if (jobAction) {
      if (jobAction.classList.contains('danger') && !confirm(jobAction.textContent + '?')) return;
      m.setv({ value: jobAction.dataset.action + ':' + jobAction.dataset.uid + ':' + Date.now() });
      return;
    }
    var btn = ev.target.closest('.object-btn');
    if (btn) {
      // losing an in-flight synthesis or fit costs real work, so those rows confirm
      // (data-confirm); a file load or a line-strength band cancels straight away.
      // dismiss-job drops a retained (failed) Processes row.
      if (btn.dataset.action === 'cancel-pending' && btn.dataset.confirm === '1' &&
          !confirm('Cancel this run?')) return;
      m.setv({ value: btn.dataset.action + ':' + btn.dataset.uid + ':' + Date.now() });
      return;
    }
    var item = ev.target.closest('.object-item');
    if (item && item.dataset.uid) {
      m.setv({ value: 'select:' + item.dataset.uid + ':' + Date.now() });
    }
  });

  // the Mask tool: toggle from the toolbar, leave via × or Esc; the panel's actions are
  // rendered dynamically, so they are delegated
  document.getElementById('mask-btn').addEventListener('click', function () { ewAction('mask-toggle'); });
  document.getElementById('range-close').addEventListener('click', function () { ewAction('range-cancel'); });
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.range-new-mask')) ewAction('range-new-mask');
    else if (ev.target.closest('.range-synth')) ewAction('range-synth');
    else if (ev.target.closest('.range-fit')) ewAction('range-fit');
    else if (ev.target.closest('.range-ew')) ewAction('range-ew');
    else if (ev.target.closest('.range-add')) ewAction('range-add');
  });

  // line-position indicators: show/hide whatever the open view attached
  document.getElementById('annot-btn').addEventListener('click', function () { ewAction('annot-toggle'); });
  document.getElementById('linelist-close').addEventListener('click', function () { ewAction('lineclose'); });

  // the line-list search: every keystroke is a query (after a short pause, so a word typed
  // at speed costs one search), Esc empties it and so returns the window's lines
  (function () {
    var box = document.getElementById('linelist-search'), timer = null;
    function send() {
      timer = null;
      var m = model('linelist-search');
      if (m) m.setv({ value: box.value + '\u0000' + Date.now() });
    }
    box.addEventListener('input', function () {
      if (timer) clearTimeout(timer);
      timer = setTimeout(send, 150);
    });
    box.addEventListener('keydown', function (ev) {
      if (ev.key === 'Escape') { box.value = ''; send(); box.blur(); ev.stopPropagation(); }
      else if (ev.key === 'Enter') { if (timer) clearTimeout(timer); send(); }
    });
  })();

  // the line list is re-rendered wholesale into its slot, so the row to jump to can only be
  // placed once it is in the DOM: the server marks it .jump-target and this positions it.
  // What to do is the server's call, carried on the row as data-scroll. "ifneeded" is a pick
  // made out on the spectrum: scroll only when the row is actually out of view. "keep" is a
  // click on a row of this table, which narrows the window and so re-cuts the table under
  // the pointer: put the row back at the height it was clicked at, so the pointer is still
  // over the same line. Every other render ("none") leaves the scroll alone.
  // scrollTop is set directly rather than via scrollIntoView, which would also scroll every
  // ancestor and shift the whole workspace.
  (function () {
    var slot = document.getElementById('linelist-table');
    if (!slot) return;
    new MutationObserver(function () {
      var row = slot.querySelector('.line-row.jump-target');
      if (!row) return;
      var mode = row.dataset.scroll;
      var box = slot.getBoundingClientRect(), r = row.getBoundingClientRect();
      if (mode === 'keep') {
        // only for the row this anchor was taken from: a jump can be clamped to a window
        // that does not hold the clicked line at all
        var anchor = linelistAnchor;
        linelistAnchor = null;
        if (!anchor || anchor.lam0 !== row.dataset.lam0
            || anchor.ion !== row.dataset.ion) return;
        // clamped by the scroll extent: a row near either end cannot always be put back
        slot.scrollTop += (r.top - box.top) - anchor.off;
        return;
      }
      if (mode !== 'ifneeded') return;
      var head = slot.querySelector('.linelist-table thead');
      // the header is sticky, so a row under it is covered, not visible
      var top = box.top + (head ? head.getBoundingClientRect().height : 0);
      if (r.top >= top && r.bottom <= box.bottom) return;
      slot.scrollTop += (r.top + r.height / 2) - (top + box.bottom) / 2;
    }).observe(slot, { childList: true, subtree: true });
  })();

  // noise filter: toggle the panel/feautrier view / clear the applied filter
  document.getElementById('filter-btn').addEventListener('click', function () { ewAction('filter'); });
  document.getElementById('filter-clear').addEventListener('click', function () { ewAction('filter-clear'); });
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.filter-clear-btn')) ewAction('filter-clear');
    else if (ev.target.closest('.filter-close')) ewAction('filter');
  });
  document.getElementById('filter-switch-input').addEventListener('change', function (ev) {
    var m = model('object-actions');
    if (m) m.setv({ value: 'filtertoggle:' + (ev.target.checked ? 'on' : 'off') + ':' + Date.now() });
  });
  // the Mask tool's equivalent-width settings (rendered into its panel)
  document.addEventListener('change', function (ev) {
    var m = model('object-actions');
    if (!m) return;
    var sel = ev.target.closest('.ew-method-select');
    if (sel) { m.setv({ value: 'profile:' + sel.value + ':' + Date.now() }); return; }
    var snr = ev.target.closest('.ew-snr-input');
    if (snr) { m.setv({ value: 'ew-snr:' + snr.value + ':' + Date.now() }); return; }
    var cont = ev.target.closest('.ew-cont-input');
    if (cont) { m.setv({ value: 'ew-cont:' + cont.value + ':' + Date.now() }); }
  });

  // continuum normalisation: toggle from the toolbar; init/clear/close from the panel
  document.getElementById('normalise-btn').addEventListener('click', function () { ewAction('normalise'); });
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.norm-init')) ewAction('norm-init');
    else if (ev.target.closest('.norm-clear')) ewAction('norm-clear');
    else if (ev.target.closest('.norm-cancel')) ewAction('norm-cancel');
    else if (ev.target.closest('.norm-close')) ewAction('normalise');
  });

  // line masks: toggle the panel from the toolbar; reveal the import picker; per-row actions
  // (select / delete) + save / clear go through the dedicated "linemask-actions" channel
  function linemaskAction(kind, arg) {
    var m = model('linemask-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.linemask-close')) { ewAction('linemask'); return; }
    // a line-mask object in the left sidebar: open its view, or drop the whole set
    var msdel = ev.target.closest('[data-maskset-remove]');
    if (msdel) { ev.stopPropagation(); linemaskAction('set-remove', msdel.dataset.masksetRemove); return; }
    var msitem = ev.target.closest('.maskset-item');
    if (msitem) { linemaskAction('set-select', msitem.dataset.maskset); return; }
    if (ev.target.closest('.lm-new')) { linemaskAction('new'); return; }
    if (ev.target.closest('.lm-synthesize')) { linemaskAction('synthesize'); return; }
    if (ev.target.closest('.lm-fit')) { linemaskAction('fit'); return; }
    if (ev.target.closest('.lm-measure-abund')) { linemaskAction('measure-abund'); return; }
    if (ev.target.closest('.lm-fit-balance')) { linemaskAction('fit-balance'); return; }
    if (ev.target.closest('.lm-revert')) {
      if (confirm('Discard every edit and reload this set from its file?')) linemaskAction('revert');
      return;
    }
    if (ev.target.closest('.lm-saveas')) {
      // the file name rides along as the action arg (space-delimited channel, so a name
      // with spaces survives — only the trailing nonce is stripped server-side)
      var nameEl = document.querySelector('.lm-name');
      linemaskAction('save-as', nameEl ? nameEl.value.trim() : '');
      return;
    }

    var check = ev.target.closest('.lm-table .check-head-btn[data-checkcol]');
    if (check) { linemaskAction('check-' + check.dataset.checkact, check.dataset.checkcol); return; }
    var del = ev.target.closest('.lm-del');
    if (del) { ev.stopPropagation(); linemaskAction('del', del.dataset.idx); return; }
    // both tables' rows carry the same data-idx, so either selects the mask
    var row = ev.target.closest('.lm-row');
    if (row && !ev.target.closest('.lm-incl')) {
      linemaskAction('select', row.dataset.idx);
    }
  });
  // the Empty-set pane: create a blank set instead of importing a file
  function newMaskSet() {
    var inp = document.getElementById('linemask-new-name');
    linemaskAction('new-set', JSON.stringify({name: (inp.value || '').trim(), element: document.getElementById('linemask-new-element').value.trim()}));
    inp.value = '';
  }
  document.getElementById('linemask-new-set').addEventListener('click', newMaskSet);
  document.getElementById('linemask-new-name').addEventListener('keydown', function (e) {
    if (e.key === 'Enter') { e.preventDefault(); newMaskSet(); }
  });

  // the Auto-generate pane: line-list picker, spectrum choice, and the run button. The
  // pane is only rendered while that mask source exists (masktool.AUTO_MASKS)
  var autoBrowse = document.getElementById('linemask-auto-browse');
  if (autoBrowse) {
    autoBrowse.addEventListener('click', function () {
      var panel = document.getElementById('linemask-auto-picker');
      panel.classList.toggle('open');
      if (panel.classList.contains('open')) {
        seedResults('linemask-auto-results');
      }
    });
    document.getElementById('linemask-auto-run').addEventListener('click', function () {
      linemaskAction('auto-run');
    });
  }

  document.addEventListener('change', function (ev) {
    var inc = ev.target.closest('.lm-incl');
    if (inc) { linemaskAction('incl', inc.dataset.idx); return; }
    var ref = ev.target.closest('.lm-auto-ref');
    if (ref) { linemaskAction('auto-ref', ref.value); }
  });

  // keyboard: a saved-workspace row activates on Enter / Space; Escape with no dialog
  // open leaves whatever tool holds the canvas
  document.addEventListener('keydown', function (ev) {
    if ((ev.key === 'Enter' || ev.key === ' ') && ev.target.matches('.session-item')) { ev.preventDefault(); ev.target.click(); }
    if (ev.key === 'Escape' && !document.querySelector('dialog[open]')) ewAction('escape');
  });

  // ---- draggable splitter for the centre bottom panel (the sidebars are fixed-width) ----
  // Layout sizes are pure client-side UI state, kept in localStorage (no server round-trips
  // per mousemove). During a drag, body.resizing sets pointer-events:none on everything so the
  // Bokeh canvas can't swallow the mousemove — the document-level handlers below own it.
  (function () {
    var LS = 'tumult.layout';
    var VARS = { 'bottom': '--bottom-h', 'linelist': '--linelist-h' };
    var saved = {};
    try { saved = JSON.parse(localStorage.getItem(LS) || '{}'); } catch (e) {}
    Object.keys(saved).forEach(function (k) {
      if (VARS[k]) document.documentElement.style.setProperty(VARS[k], saved[k]);
    });
    function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }
    function set(kind, value) {
      document.documentElement.style.setProperty(VARS[kind], value);
      saved[kind] = value;
    }
    var drag = null;
    document.addEventListener('mousedown', function (e) {
      var sp = e.target.closest('.splitter');
      if (!sp) return;
      e.preventDefault();
      drag = sp.dataset.resize;
      document.body.classList.add('resizing');
    });
    document.addEventListener('mousemove', function (e) {
      if (!drag) return;
      var m = document.getElementById('main').getBoundingClientRect();
      // both bottom panels are measured from the bottom of the centre column; the line list
      // is the last child, so it is the one the raw distance describes when both are open
      if (drag === 'bottom' || drag === 'linelist') {
        set(drag, clamp((m.bottom - e.clientY) / m.height * 100, 15, 85) + '%');
      }
    });
    document.addEventListener('mouseup', function () {
      if (!drag) return;
      drag = null;
      document.body.classList.remove('resizing');
      try { localStorage.setItem(LS, JSON.stringify(saved)); } catch (e) {}
      window.dispatchEvent(new Event('resize'));  // nudge Bokeh figures to reflow to the new size
    });
  })();
})();
