// TumultStudio preferences: the Preferences modal — its category rail, the run-cache
// page, the atmosphere-grid cards with their shared folder picker, the line-list
// downloads and the reference-star table. Every row is rendered by the server (preferences/), so the
// handlers are delegated. Binds after bridge.js and workspace.js.
(function () {
  var TS = window.TumultStudio;
  var model = TS.model, wireModal = TS.wireModal, wireNav = TS.wireNav, tailPaths = TS.tailPaths;

  wireModal('open-settings', 'close-settings', 'settings-modal');

  // cache page actions go through the dedicated "cache-actions" channel
  function cacheAction(kind, arg) {
    var m = model('cache-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }

  wireNav('#settings-modal', function (page) {
    if (page === 'cache') cacheAction('refresh');  // re-scan disk when opened
    if (page === 'linelists') linelistAction('list');  // list the share once, when opened
  });

  // clear-all + per-row delete (rows are rendered dynamically into #synth-cache)
  document.getElementById('cache-clear').addEventListener('click', function () {
    if (confirm('Delete all cached synthesis runs?')) cacheAction('clear');
  });
  document.addEventListener('click', function (ev) {
    var del = ev.target.closest('.cache-del');
    if (del && confirm('Delete this cached run from disk?')) cacheAction('del', del.dataset.path);
  });

  // atmosphere grids (Preferences → Atmosphere grids)
  function gridAction(kind, arg) {
    var m = model('grid-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }
  // edit a field of the storage settings or a custom grid card -> "grid-fields-input"
  document.addEventListener('change', function (ev) {
    var el = ev.target.closest('[data-grid-key]');
    if (!el) return;
    var m = model('grid-fields-input');
    if (m) m.setv({ value: el.dataset.gridId + '|' + el.dataset.gridKey + '=' + el.value + ' ' + Date.now() });
  });
  // pick the 1D LTE reference grid (for the abundance-correction tool) -> "grid-actions"
  document.addEventListener('change', function (ev) {
    if (ev.target.closest('.grid-ref-select')) gridAction('setref', ev.target.value);
  });
  // the folder row of a grid card is the folder picker: focusing its field hands the row to
  // the one shared picker (#grid-picker, moved into the row in place of the plain field),
  // whose listing drops down under it and leads with a "use this folder" row. Leaving the
  // picker without choosing gives the row its plain field back. The card list is
  // re-rendered from the server on every edit, which detaches the picker, so it is put back
  // after each render (the node survives; Bokeh's view stays bound to it).
  var gridPicker = document.getElementById('grid-picker'), gridList = document.getElementById('grid-list');
  function gridPickerInput() {
    var found = null;
    (function walk(root) {
      root.querySelectorAll('*').forEach(function (n) {
        if (found || !n.shadowRoot) return;
        found = n.shadowRoot.querySelector('input.bk-input') || null;
        if (!found) walk(n.shadowRoot);
      });
    })(gridPicker);
    return found;
  }
  function placeGridPicker() {
    var row = gridPicker.dataset.gridId
      && gridList.querySelector('.grid-folder-row[data-grid-id="' + gridPicker.dataset.gridId + '"]');
    gridList.querySelectorAll('.grid-folder-row.editing').forEach(function (r) {
      if (r !== row) r.classList.remove('editing');
    });
    if (row) { row.classList.add('editing'); row.appendChild(gridPicker); }
    else { gridPicker.classList.remove('open'); gridList.after(gridPicker); }
  }
  function leaveGridPicker() {
    gridPicker.dataset.gridId = '';
    placeGridPicker();
  }
  new MutationObserver(function () {
    if (!gridPicker.isConnected) placeGridPicker();
    tailPaths(gridList);
  }).observe(gridList, { childList: true });
  // the server closes the picker once a folder is adopted
  new MutationObserver(function () {
    if (gridPicker.dataset.gridId && !gridPicker.classList.contains('open')) leaveGridPicker();
  }).observe(gridPicker, { attributes: true, attributeFilter: ['class'] });
  document.addEventListener('focusin', function (ev) {
    var field = ev.target.closest && ev.target.closest('input.grid-folder');
    if (!field) return;
    gridPicker.dataset.gridId = field.dataset.gridId;
    placeGridPicker();
    gridPicker.classList.add('open');
    gridAction('browse', field.dataset.gridId);   // the picker lists the row's folder
    var inp = gridPickerInput();
    if (inp) inp.focus();
  });
  document.addEventListener('focusout', function (ev) {
    if (!gridPicker.dataset.gridId || !ev.composedPath().includes(gridPicker)) return;
    // focus moving within the picker (its list keeps focus via mousedown) is not a leave
    setTimeout(function () {
      var active = document.activeElement;
      if (!gridPicker.contains(active)) leaveGridPicker();
    }, 0);
  });
  // the grid rows, model rows and custom cards (rendered dynamically into #grid-list):
  // "<kind> <grid id> [<model>]"; deleting files asks first
  document.addEventListener('click', function (ev) {
    var sortth = ev.target.closest('#grid-list .grid-models .line-sort');
    if (sortth) { gridAction('sort', sortth.dataset.sortcol); return; }
    var btn = ev.target.closest('#grid-list button');
    if (!btn || btn.disabled) return;
    var gid = btn.dataset.gridId || '', model_ = btn.dataset.model || '';
    if (btn.classList.contains('grid-add')) gridAction('add');
    else if (btn.classList.contains('grid-del')) gridAction('del', gid);
    else if (btn.classList.contains('grid-back')) gridAction('back');
    else if (btn.classList.contains('grid-edit')) gridAction('edit', gid);
    else if (btn.classList.contains('grid-fetch-all')) gridAction('fetch-all', gid);
    else if (btn.classList.contains('grid-fetch-model')) gridAction('fetch-model', gid + ' ' + model_);
    else if (btn.classList.contains('grid-del-catalog')) {
      if (confirm('Delete every downloaded file of ' + btn.dataset.name + '?')) gridAction('del-catalog', gid);
    } else if (btn.classList.contains('grid-del-model')) {
      if (confirm('Delete the downloaded model ' + btn.dataset.name + '?')) gridAction('del-model', gid + ' ' + model_);
    }
  });

  // line lists (Preferences → Line lists): the head row and the file rows (rendered
  // dynamically into #linelist-list) -> "linelist-actions" as "<kind> [<file>]"; deleting
  // files asks first
  function linelistAction(kind, arg) {
    var m = model('linelist-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }
  document.addEventListener('click', function (ev) {
    var sortth = ev.target.closest('#linelist-list .line-sort');
    if (sortth) { linelistAction('sort', sortth.dataset.sortcol); return; }
    var btn = ev.target.closest('#linelist-list button');
    if (!btn || btn.disabled) return;
    var name = btn.dataset.name || '';
    if (btn.classList.contains('linelist-relist')) linelistAction('list');
    else if (btn.classList.contains('linelist-fetch-all')) linelistAction('fetch-all');
    else if (btn.classList.contains('linelist-fetch')) linelistAction('fetch', name);
    else if (btn.classList.contains('linelist-del-all')) {
      if (confirm('Delete every downloaded master line list?')) linelistAction('del-all');
    } else if (btn.classList.contains('linelist-del')) {
      if (confirm('Delete the downloaded line list ' + name + '?')) linelistAction('del', name);
    }
  });

  // reference stars (Preferences → Reference stars): edits of a row's inputs report
  // through "refstars-fields-input", the buttons through "refstars-actions"
  document.addEventListener('change', function (ev) {
    var el = ev.target.closest('[data-star-key][data-star-field]');
    if (!el) return;
    var m = model('refstars-fields-input');
    if (m) m.setv({ value: el.dataset.starKey + '|' + el.dataset.starField + '=' + el.value + ' ' + Date.now() });
  });
  document.addEventListener('click', function (ev) {
    var btn = ev.target.closest('#refstars-list button');
    if (!btn) return;
    var m = model('refstars-actions');
    if (!m) return;
    var act = function (kind, arg) { m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() }); };
    if (btn.classList.contains('refstar-add')) act('add');
    else if (btn.classList.contains('refstar-del')) act('del', btn.dataset.starKey);
    else if (btn.classList.contains('refstar-reset')) {
      if (confirm('Replace the reference stars by the default list?')) act('reset');
    }
  });
})();
