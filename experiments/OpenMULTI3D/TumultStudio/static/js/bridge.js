// TumultStudio bridge: the one place the page reaches the server's bridge models, and
// the behaviour of the controls every feature builds on (ui/): dialogs, the category
// rail, drop-down menus, file pickers, tabs, collapsed sections, the submit-button busy
// state. Loaded before the feature scripts, which read these helpers off
// window.TumultStudio and bind their own handlers after the markup they act on exists.
(function () {
  // a named bridge model of the document (see htmlbridge): null until Bokeh has it
  function model(name) { return Bokeh.documents[0].get_model_by_name(name); }

  // a file-picker box holding a path longer than itself shows the path's tail
  // (…/Marcs/model), which only takes once the box is laid out: on every modal open and
  // whenever a box is left. The boxes are Bokeh widgets, so the walk descends shadow roots.
  function tailPaths(root) {
    root.querySelectorAll('input.grid-folder').forEach(function (i) { i.scrollLeft = i.scrollWidth; });
    root.querySelectorAll('*').forEach(function (n) {
      if (!n.shadowRoot) return;
      n.shadowRoot.querySelectorAll('input.bk-input[type="text"]').forEach(function (i) {
        i.scrollLeft = i.scrollWidth;
      });
      tailPaths(n.shadowRoot);
    });
  }
  document.addEventListener('focusout', function (ev) {
    var el = ev.composedPath()[0];
    if (el && el.tagName === 'INPUT'
        && (el.classList.contains('bk-input') || el.classList.contains('grid-folder'))) {
      setTimeout(function () { el.scrollLeft = el.scrollWidth; }, 0);
    }
  }, true);
  // a box scrolls once it is laid out, so watch what reveals one: a dialog opening (from
  // the page or a server bridge alike), an inline picker gaining its open class, a
  // settings page or tab panel becoming the active one
  new MutationObserver(function (muts) {
    muts.forEach(function (m) {
      var t = m.target;
      if ((t.tagName === 'DIALOG' && t.open)
          || t.classList.contains('open') || t.classList.contains('active')) {
        setTimeout(function () { tailPaths(t); }, 0);
      }
    });
  }).observe(document.body, { attributes: true, subtree: true, attributeFilter: ['open', 'class'] });

  document.addEventListener('click', function(ev) {
    if (ev.target.closest('.modal-cancel')) ev.target.closest('dialog').close();
  });
  // collapsed advanced sections (controls.Form groups): a form re-renders on every server
  // edit, so remember each section's open state and put it back on the rebuilt markup
  var disclosures = new Map();
  document.addEventListener('toggle', function(ev) {
    if (ev.target.matches('details[data-disclosure]')) disclosures.set(ev.target.dataset.disclosure, ev.target.open);
  }, true);
  new MutationObserver(function (muts) {
    muts.forEach(function (m) {
      m.addedNodes.forEach(function (node) {
        if (node.nodeType !== 1) return;
        var found = node.matches('details[data-disclosure]') ? [node]
                  : Array.from(node.querySelectorAll('details[data-disclosure]'));
        found.forEach(function (d) { if (disclosures.get(d.dataset.disclosure)) d.open = true; });
      });
    });
  }).observe(document.body, { childList: true, subtree: true });
  // the footer's status message is mirrored into the status slot of every open dialog
  var statusObserver = new MutationObserver(function() {
    var status = document.querySelector('#footer-status');
    document.querySelectorAll('dialog[open]').forEach(function(dialog) {
      var slot = dialog.querySelector('.modal-status');
      if (slot && status) slot.textContent = status.textContent;
    });
  });
  statusObserver.observe(document.getElementById('footer-status'), {childList: true, subtree: true});
  document.addEventListener('invalid', function(ev) {
    var section = ev.target.closest('details'); if (section) section.open = true;
  }, true);
  // wire a <dialog>: `openId` may be null for a modal opened from a menu instead of a
  // dedicated button — the returned function opens it either way.
  function wireModal(openId, closeId, modalId, onOpen) {
    var m = document.getElementById(modalId);
    var opener;
    var status = document.createElement('div'); status.className = 'modal-status hint';
    status.setAttribute('role', 'status'); m.appendChild(status);
    m.addEventListener('close', function() {
      var channel = model('dialog-actions');
      if (channel) channel.value = modalId + ' ' + Date.now();
      if (opener && opener.isConnected) opener.focus();
      // a job-launching submit button left busy (see busySubmit) is only ever inside a
      // static footer that isn't re-rendered on reopen — reset it here, once, whatever
      // closed the dialog (a successful submit, Cancel, Esc, or the backdrop)
      resetBusy(m);
    });
    function open() { opener = document.activeElement; m.showModal(); if (onOpen) onOpen(); }
    var o = openId && document.getElementById(openId);
    if (o) o.addEventListener('click', open);
    var c = document.getElementById(closeId);
    if (c) c.addEventListener('click', function () { m.close(); });
    m.addEventListener('click', function (e) { if (e.target === m) m.close(); });
    return open;
  }

  // pre-render a picker's current results into its container, which names its results
  // bridge as data-results: the html bridge's js_on_change only fires on change, so a
  // freshly-opened modal would start empty
  function seedResults(containerId) {
    var el = document.getElementById(containerId);
    var rb = el && el.dataset.results && model(el.dataset.results);
    if (rb && el) el.innerHTML = rb.value;
  }

  // a left nav rail switching .settings-page panes inside one container (the Preferences
  // categories and the New line mask sources both use it)
  function wireNav(rootSel, onSwitch) {
    var root = document.querySelector(rootSel);
    if (!root) return;
    root.querySelectorAll('.settings-nav-item').forEach(function (tab) {
      tab.addEventListener('click', function () {
        var page = tab.dataset.page;
        root.querySelectorAll('.settings-nav-item').forEach(function (t) {
          t.classList.toggle('active', t === tab);
        });
        root.querySelectorAll('.settings-page').forEach(function (p) {
          p.classList.toggle('active', p.dataset.page === page);
        });
        if (onSwitch) onSwitch(page);
      });
    });
  }

  // inline atmosphere file pickers (synth / abund / corr, the preference paths): the
  // "Type a path…" field is always shown while a concrete file is in play; its result
  // list drops down while the field is focused and — crucially — stays open while
  // folders in the list are clicked (mousedown inside the list is prevented, so the field
  // keeps focus and navigates in place). Clicking away closes it; picking a file
  // collapses it from the server.
  // The list itself is a top-layer popover (see base.css) so it can overhang the column
  // and the modal: its open state stays the `open` + `results-open` class pair (set here
  // and pulsed from the server), and syncDropdown keeps the popover in step with them.
  function placeDropdown(pick, list) {
    var r = pick.getBoundingClientRect();
    list.style.left = r.left + 'px';
    list.style.top = (r.bottom + 2) + 'px';
    list.style.minWidth = r.width + 'px';
    list.style.maxWidth = (window.innerWidth - r.left - 8) + 'px';
    list.style.maxHeight = Math.max(80, Math.min(200, window.innerHeight - r.bottom - 12)) + 'px';
  }
  function syncDropdown(pick) {
    var list = pick.querySelector('.file-results');
    if (!list || !list.showPopover) return;
    var want = pick.classList.contains('open') && pick.classList.contains('results-open')
               && pick.getClientRects().length > 0;
    var shown = list.matches(':popover-open');
    if (want) { placeDropdown(pick, list); if (!shown) list.showPopover(); }
    else if (shown) list.hidePopover();
  }
  function syncDropdowns() {
    document.querySelectorAll('.inline-picker-auto').forEach(syncDropdown);
  }
  document.addEventListener('focusin', function (ev) {
    var pick = ev.target.closest('.inline-picker-auto');
    if (!pick) return;
    pick.classList.add('results-open');
    var el = pick.querySelector('.file-results'), nm = pick.dataset.results;
    if (el && nm) { var rb = model(nm); if (rb) el.innerHTML = rb.value; }
    syncDropdown(pick);
  });
  document.addEventListener('focusout', function (ev) {
    var pick = ev.target.closest('.inline-picker-auto');
    if (pick) { pick.classList.remove('results-open'); syncDropdown(pick); }
  });
  // the server pulses the classes; a modal closing or a page scrolling moves the field
  new MutationObserver(function (muts) {
    muts.forEach(function (m) {
      if (m.target.classList.contains('inline-picker-auto')) syncDropdown(m.target);
    });
  }).observe(document.body, { attributes: true, subtree: true, attributeFilter: ['class'] });
  document.addEventListener('scroll', syncDropdowns, true);
  window.addEventListener('resize', syncDropdowns);
  document.addEventListener('close', function (ev) {
    if (ev.target.tagName === 'DIALOG') syncDropdowns();
  }, true);
  document.addEventListener('mousedown', function (ev) {
    var pick = ev.target.closest('.inline-picker-auto');
    if (!pick) return;
    if (ev.target.closest('.file-results')) ev.preventDefault();  // keep focus → list stays open
    else pick.classList.add('results-open');                      // clicking the field (re)opens
  });

  // a drop-down menu under a toolbar button: open/close, arrow-key and Escape handling,
  // and `action(item)` for the item picked
  function wireMenu(buttonId, menuId, action, onOpen) {
    var button = document.getElementById(buttonId), menu = document.getElementById(menuId);
    // read on use: a server-rendered menu (Edit) changes its items
    function enabledItems() {
      return Array.from(menu.querySelectorAll('[role="menuitem"]:not([disabled])'));
    }
    function close(focus) {
      menu.classList.remove('open'); button.setAttribute('aria-expanded', 'false');
      if (focus) button.focus();
    }
    button.addEventListener('click', function(ev) {
      ev.stopPropagation();
      var open = !menu.classList.contains('open');
      document.querySelectorAll('.menu.open').forEach(function(m) { m.classList.remove('open'); });
      document.querySelectorAll('.menu-btn').forEach(function(b) { b.setAttribute('aria-expanded', 'false'); });
      menu.classList.toggle('open', open); button.setAttribute('aria-expanded', String(open));
      if (open && onOpen) onOpen();
      var items = enabledItems();
      if (open) (items[0] || menu).focus();
    });
    button.addEventListener('keydown', function(ev) { if (ev.key === 'ArrowDown') { ev.preventDefault(); button.click(); } });
    menu.addEventListener('keydown', function(ev) {
      var items = enabledItems(), i = items.indexOf(document.activeElement);
      if (ev.key === 'Escape') { ev.preventDefault(); ev.stopPropagation(); close(true); }
      else if (ev.key === 'ArrowDown' || ev.key === 'ArrowUp') {
        if (!items.length) return;
        ev.preventDefault(); items[(i + (ev.key === 'ArrowDown' ? 1 : items.length - 1)) % items.length].focus();
      } else if (ev.key === 'Tab') close(false);
    });
    menu.addEventListener('click', function(ev) {
      var item = ev.target.closest('.menu-item'); if (!item || item.disabled) return;
      close(true); action(item);
    });
    document.addEventListener('click', function() { close(false); });
  }

  // file pickers (filepicker.FilePicker): a row names its picker bridge and path. In the
  // Open dialogs a click selects and Open / double-click / Enter picks; everywhere else
  // (inline pickers) a click picks straight away, as does a folder row anywhere.
  function pickFile(li) {
    var m = model(li.dataset.picker);
    if (m) m.value = li.dataset.path + '\u0000' + Date.now();
  }
  function isOpenPicker(li) {
    var dialog = li.closest('dialog');
    return dialog && ['file-modal', 'atom-modal', 'linemask-modal'].includes(dialog.id);
  }
  document.addEventListener('click', function(ev) {
    var li = ev.target.closest('.file-item');
    if (li) {
      if (!isOpenPicker(li) || li.dataset.dir === '1') { pickFile(li); return; }
      var dialog = li.closest('dialog');
      dialog.querySelectorAll('.file-item.selected').forEach(function(r) { r.classList.remove('selected'); });
      li.classList.add('selected');
      var button = dialog.querySelector('.picker-open'); if (button) button.disabled = false;
    }
    var open = ev.target.closest('.picker-open');
    if (open) {
      var selected = open.closest('dialog').querySelector('.file-item.selected');
      if (selected) pickFile(selected);
    }
  });
  document.addEventListener('dblclick', function(ev) {
    var li = ev.target.closest('.file-item');
    if (li && isOpenPicker(li) && li.dataset.dir !== '1') pickFile(li);
  });
  document.addEventListener('keydown', function(ev) {
    var li = ev.target.closest('.file-item'); if (!li) return;
    if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); pickFile(li); }
    if (ev.key === 'ArrowDown' || ev.key === 'ArrowUp') {
      ev.preventDefault(); var rows = Array.from(li.parentNode.children), i = rows.indexOf(li);
      rows[(i + (ev.key === 'ArrowDown' ? 1 : rows.length - 1)) % rows.length].focus();
    }
  });

  // generic tab switching, scoped per group (a .tab-bar and its adjacent .tab-content)
  document.querySelectorAll('.tab-bar').forEach(function (bar) {
    // a panel may have more than one .tab-content sibling — switch the active panel in
    // every one of them together
    var contents = bar.parentElement.querySelectorAll(':scope > .tab-content');
    bar.querySelectorAll('.tab').forEach(function (tab) {
      tab.addEventListener('click', function () {
        bar.querySelectorAll('.tab').forEach(function (t) { t.classList.toggle('active', t === tab); });
        contents.forEach(function (content) {
          content.querySelectorAll('.tab-panel').forEach(function (p) {
            p.classList.toggle('active', p.dataset.panel === tab.dataset.panel);
          });
        });
      });
    });
  });

  // named actions on the "object-actions" channel (tool toggles, panel buttons, …)
  function ewAction(name) {
    var m = model('object-actions');
    if (m) m.setv({ value: name + ':-:' + Date.now() });
  }

  // a job-launching modal's submit button: the daemon round trip a click starts can take
  // a moment, so give feedback at once (disabled, with an animated ellipsis after the
  // label — the .busy rule in base.css) instead of leaving the button looking inert until
  // the server gets around to closing the dialog. Returns false (and sends nothing) for a
  // click that arrives while the button is already busy — the client-side half of the
  // double-submit guard; the synth modal's own `_pending` flag (synthmodal.py) holds the
  // server-side half for its action. A rejected submit (the dialog stays open) resets the
  // button through the dialog bridge's "reject" state (see htmlbridge.dialog_bridge)
  // rather than here, since that is a server decision.
  function busySubmit(btn) {
    if (!btn || btn.disabled) return false;
    btn.disabled = true;
    btn.classList.add('busy');
    return true;
  }
  function resetBusy(scope) {
    (scope || document).querySelectorAll('button.busy').forEach(function (btn) {
      btn.classList.remove('busy');
      btn.disabled = false;
    });
  }

  window.TumultStudio = {
    model: model, tailPaths: tailPaths, wireModal: wireModal, seedResults: seedResults,
    wireNav: wireNav, wireMenu: wireMenu, ewAction: ewAction,
    busySubmit: busySubmit, resetBusy: resetBusy,
  };
})();
