// TumultStudio analysis: the abundance, NLTE-correction, model-atom and
// ionisation-balance tools — each tool's run and close wiring on its own
// "<prefix>-actions" channel, the abundance modal's model-atom selector, and the
// table <-> scatter hover link of the diagnostic panels. The interactions the tools
// share (step / select / remove / sort / row jump / include toggles) are bound once by
// resultview's delegated script, emitted by the page after these files and keyed by
// the data-tool scope of the markup.
(function () {
  var TS = window.TumultStudio;
  var model = TS.model, wireModal = TS.wireModal, seedResults = TS.seedResults;
  var ewAction = TS.ewAction, busySubmit = TS.busySubmit;

  // the abundance modal's model-atom selector (NLTE): reveal the picker, delete chips
  document.getElementById('browse-abund-atom').addEventListener('click', function() {
    document.getElementById('abund-atom-picker').classList.toggle('open');
    seedResults('abund-atom-results');
  });
  document.addEventListener('click', function(ev) {
    var del = ev.target.closest('.abund-atom-del'); if (!del) return;
    var m = model('abund-atom-actions');
    if (m) m.value = 'del ' + del.dataset.path + ' ' + Date.now();
  });

  // Table <-> scatter hover link. A row and a point carry the same line index (the row as
  // data-line, the point as the source's `idx` column), so the marker source published as
  // "<prefix>_hover" can be moved onto the hovered row's point without a server round trip.
  // The other direction is the scatter's own HoverTool callback (see diagnostics.HOVER_JS).
  function markPoint(prefix, line) {
    var doc = Bokeh.documents[0];
    var hov = doc.get_model_by_name(prefix + '_hover');
    if (!hov) return;
    var empty = { x_evi: [], x_rew: [], a: [] };
    var names = [prefix + '_points', prefix + '_points2'];
    for (var n = 0; n < names.length; n++) {
      var src = doc.get_model_by_name(names[n]);
      if (!src || !src.data.idx) continue;
      var k = Array.prototype.indexOf.call(src.data.idx, line);
      if (k >= 0) {
        hov.data = { x_evi: [src.data.x_evi[k]], x_rew: [src.data.x_rew[k]],
                     a: [src.data.a[k]] };
        return;
      }
    }
    hov.data = empty;
  }
  function clearMarks(prefix) {
    markPoint(prefix, -1);
    document.querySelectorAll('.' + prefix + '-row.row-hover')
            .forEach(function (r) { r.classList.remove('row-hover'); });
  }
  ['abund', 'corr', 'ionbal'].forEach(function (prefix) {
    var table = document.getElementById(prefix + '-table');
    var pane = document.getElementById(prefix + '-scatter');
    if (table) {
      table.addEventListener('mouseover', function (ev) {
        var row = ev.target.closest('.' + prefix + '-row');
        if (row) markPoint(prefix, +row.dataset.line);
      });
      table.addEventListener('mouseleave', function () { markPoint(prefix, -1); });
    }
    // the hover callback fires while the cursor is over the plot; leaving it never reports
    // an empty hit, so the row mark is dropped here
    if (pane) pane.addEventListener('mouseleave', function () { clearMarks(prefix); });
  });

  // each tool's own run and close wiring, on the same "<prefix>-actions" channel the
  // shared interactions use
  function analysisAction(prefix, kind, arg) {
    var m = model(prefix + '-actions');
    if (m) m.setv({ value: kind + ':' + (arg || '') + ':' + Date.now() });
  }

  // abundance tool: the panel's run buttons, the settings modal and the reference select
  wireModal(null, 'close-abund', 'abund-modal');
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.abund-close')) { ewAction('abundance'); return; }
    if (ev.target.closest('.abund-settings')) { analysisAction('abund', 'settings'); return; }
    if (ev.target.closest('.abund-run-lines')) { analysisAction('abund', 'run-lines'); return; }
    if (ev.target.closest('.abund-run-synth')) { analysisAction('abund', 'run-synth'); return; }
    if (ev.target.closest('.abund-run-uncert')) { analysisAction('abund', 'run-uncert'); }
  });
  document.addEventListener('change', function (ev) {
    var ref = ev.target.closest('.abund-ref-select');
    if (ref) analysisAction('abund', 'ref', ref.value);
  });

  // abundance-correction tool: the modal's Compute button
  wireModal(null, 'close-corr', 'corr-modal');
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.abundcorr-close')) { ewAction('abundcorr'); return; }
    var run = ev.target.closest('.corr-run');
    if (run && busySubmit(run)) analysisAction('corr', 'run');
  });

  // model atoms: the left-sidebar objects (select / remove), the panel's corrections action,
  // the view's tab switch, header sorts and the per-table search fields. Space-delimited so a
  // query keeps its spaces; only the trailing nonce is stripped server-side.
  function atomAction(kind, arg) {
    var m = model('atomobj-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.atom-close')) { ewAction('atom'); return; }
    if (ev.target.closest('.atom-corr')) { atomAction('corr'); return; }
    var rem = ev.target.closest('[data-atom-remove]');
    if (rem) {
      ev.stopPropagation();
      // the server asks first when correction results go with it (removal.py)
      atomAction('remove', rem.dataset.atomRemove);
      return;
    }
    var aitem = ev.target.closest('.atom-item');
    if (aitem) { atomAction('select', aitem.dataset.atom); return; }
    var tab = ev.target.closest('.atom-tab-btn');
    if (tab) { atomAction('tab', tab.dataset.tab); return; }
    var asort = ev.target.closest('.atom-table .line-sort');
    if (asort) {
      var wrap = asort.closest('.atom-table-wrap');
      atomAction('sort', (wrap.id === 'atom-trans-table' ? 'trans' : 'levels') + ' ' + asort.dataset.sortcol);
    }
  });
  document.querySelectorAll('.atom-search').forEach(function (label) {
    var box = label.querySelector('input'), timer = null;
    function send() { timer = null; atomAction('search', label.dataset.tab + ' ' + box.value); }
    box.addEventListener('input', function () {
      if (timer) clearTimeout(timer);
      timer = setTimeout(send, 150);
    });
    box.addEventListener('keydown', function (ev) {
      if (ev.key === 'Escape') { box.value = ''; send(); box.blur(); ev.stopPropagation(); }
      else if (ev.key === 'Enter') { if (timer) clearTimeout(timer); send(); }
    });
  });

  // ionisation balance: the modal's run button and the panel's refit. The tool has no
  // toolbar button — a measurement object in the left sidebar is what opens it.
  wireModal(null, 'close-ionbal', 'ionbal-modal');
  document.getElementById('ionbal-run').addEventListener('click', function () {
    if (busySubmit(this)) analysisAction('ionbal', 'run');
  });
  document.addEventListener('click', function (ev) {
    if (ev.target.closest('.ionbal-close')) { ewAction('ionbal'); return; }
    if (ev.target.closest('.ionbal-refit')) { analysisAction('ionbal', 'refit'); }
  });
})();
