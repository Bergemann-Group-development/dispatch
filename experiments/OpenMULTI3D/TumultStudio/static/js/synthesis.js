// TumultStudio synthesis: the synthesise / fit modal — its model-atom selector, the Fit
// column's abundance rows and the submit button. The atmosphere selects of every
// calculation setup are driven by synthsetup's own delegated script (select[data-atmos]
// -> "<prefix>-grid-actions"), emitted by the page after these files.
(function () {
  var TS = window.TumultStudio;
  var model = TS.model, wireModal = TS.wireModal, seedResults = TS.seedResults;
  var ewAction = TS.ewAction, busySubmit = TS.busySubmit;

  // model-atom selector (NLTE): reveal the inline picker, and delete chips
  document.getElementById('browse-atom').addEventListener('click', function () {
    var panel = document.getElementById('atom-picker');
    panel.classList.toggle('open');
    if (panel.classList.contains('open')) {
      seedResults('atom-results');  // js_on_change only fires on change
    }
  });
  document.addEventListener('click', function (ev) {
    var del = ev.target.closest('.atom-del');
    if (!del) return;
    var m = model('atom-actions');
    if (m) m.setv({ value: 'del ' + del.dataset.path + ' ' + Date.now() });
  });

  // the modal's Fit section: per-row abundance removal goes through the "fit-actions" channel
  function fitAction(kind, arg) {
    var m = model('fit-actions');
    if (m) m.setv({ value: kind + ' ' + (arg || '') + ' ' + Date.now() });
  }
  document.addEventListener('click', function (ev) {
    var del = ev.target.closest('.fit-abund-del');
    if (del) fitAction('del', del.dataset.el);
  });

  // the synthesis modal: close via × / backdrop, submit via the footer button
  wireModal(null, 'close-synth', 'synth-modal');
  document.addEventListener('click', function (ev) {
    var btn = ev.target.closest('.synth-submit');
    if (btn && busySubmit(btn)) ewAction('synth-submit');
  });
})();
