"""A tiny server->client bridge: a hidden widget whose value is pushed as
innerHTML into a light-DOM element (so it can be styled with ordinary CSS and
support event delegation, unlike Bokeh's shadow-DOM widgets)."""
import html
import itertools

from bokeh.models import CustomJS, TextInput

_nonce = itertools.count()


def html_bridge(container_id, name=None):
    """Return a hidden TextInput; setting its `.value` replaces the innerHTML of
    the element with id `container_id`."""
    bridge = TextInput(value="", visible=False)
    if name:
        bridge.name = name
    bridge.js_on_change("value", CustomJS(args=dict(cid=container_id), code="""
        const el = document.getElementById(cid);
        if (el) el.innerHTML = cb_obj.value;
    """))
    return bridge


def html_bridge_all(selector, name=None):
    """Like :func:`html_bridge`, but pushes the value into the innerHTML of *every*
    element matching the CSS `selector` (e.g. the same content shown in a modal and
    inline at once)."""
    bridge = TextInput(value="", visible=False)
    if name:
        bridge.name = name
    bridge.js_on_change("value", CustomJS(args=dict(sel=selector), code="""
        document.querySelectorAll(sel).forEach(function (el) {
            el.innerHTML = cb_obj.value;
        });
    """))
    return bridge


def check_bridge(input_id):
    """Return a hidden TextInput; setting `.value` to "on…"/"off…" sets the `checked`
    state of the checkbox with id `input_id`."""
    bridge = TextInput(value="off", visible=False)
    bridge.js_on_change("value", CustomJS(args=dict(cid=input_id), code="""
        const el = document.getElementById(cid);
        if (el) el.checked = cb_obj.value.split(':')[0] === 'on';
    """))
    return bridge


def dialog_bridge(dialog_id, name=None):
    """Return a hidden TextInput; setting `.value` to "open…"/"close…" shows or closes the
    ``<dialog>`` with id `dialog_id` — so the server can dismiss a modal once the work it
    was opened for (e.g. importing the picked file) has actually succeeded. A third state,
    "reject…", leaves the dialog open (server-side validation refused the submit) but
    still resets a submit button a click left disabled and spinning — the page's own
    ``resetBusy(el)``, called if it exists, so a plain page never has to define it."""
    bridge = TextInput(value="", visible=False)
    if name:
        bridge.name = name
    bridge.js_on_change("value", CustomJS(args=dict(did=dialog_id), code="""
        const el = document.getElementById(did);
        if (!el) return;
        const want = cb_obj.value.split(':')[0];
        if (want === 'open' && !el.open) el.showModal();
        else if (want === 'close' && el.open) el.close();
        else if (want === 'reject' && typeof resetBusy === 'function') resetBusy(el);
    """))
    return bridge


def value_bridge(input_id):
    """Return a hidden TextInput; setting `.value` to ``text + NUL + nonce`` writes ``text``
    into the ``.value`` of the ``<input>`` with id `input_id` — for a field the page edits
    freely but the server occasionally resets (a search box it has stopped searching for).
    NUL separates text and nonce because the text itself may hold anything else."""
    bridge = TextInput(value="", visible=False)
    bridge.js_on_change("value", CustomJS(args=dict(cid=input_id), code="""
        const el = document.getElementById(cid);
        if (el) el.value = cb_obj.value.split('\\0')[0];
    """))
    return bridge


def confirm_bridge(channel, name=None):
    """A hidden TextInput; setting ``.value`` to ``json([message, reply]) + NUL + nonce``
    shows the browser's confirmation with ``message`` and, when accepted, writes
    ``reply`` (plus a nonce) to the action channel named ``channel``."""
    bridge = TextInput(value="", visible=False)
    if name:
        bridge.name = name
    bridge.js_on_change("value", CustomJS(args=dict(channel=channel), code="""
        const cut = cb_obj.value.indexOf('\\0');
        if (cut < 0) return;
        const [message, reply] = JSON.parse(cb_obj.value.slice(0, cut));
        if (!window.confirm(message)) return;
        const m = Bokeh.documents[0].get_model_by_name(channel);
        if (m) m.setv({ value: reply + ' ' + Date.now() });
    """))
    return bridge


def class_bridge(selector, klass):
    """Return a hidden TextInput; setting `.value` to "on…"/"off…" toggles `klass`
    on the element matched by `selector`."""
    bridge = TextInput(value="off", visible=False)
    bridge.js_on_change("value", CustomJS(args=dict(sel=selector, k=klass), code="""
        const el = document.querySelector(sel);
        if (el) el.classList.toggle(k, cb_obj.value.split(':')[0] === 'on');
    """))
    return bridge


class BridgeSet:
    """Creates the hidden server->client bridges and remembers each one.

    ``main.build`` makes ~30 of these, and every one must be added to a document root
    or the control it drives silently stops updating. Collecting them as they are made
    (``set.models``) removes the long, easy-to-forget hand-written list in the root
    ``column(...)``. Each factory mirrors the module-level function of the same purpose.
    """

    def __init__(self):
        self.models = []

    def _keep(self, bridge):
        self.models.append(bridge)
        return bridge

    def html(self, container_id, name=None):
        return self._keep(html_bridge(container_id, name=name))

    def html_all(self, selector, name=None):
        return self._keep(html_bridge_all(selector, name=name))

    def cls(self, selector, klass):
        return self._keep(class_bridge(selector, klass))

    def dialog(self, dialog_id, name=None):
        return self._keep(dialog_bridge(dialog_id, name=name))

    def check(self, input_id):
        return self._keep(check_bridge(input_id))

    def value(self, input_id):
        return self._keep(value_bridge(input_id))

    def confirm(self, channel, name=None):
        return self._keep(confirm_bridge(channel, name=name))

    def action(self, name):
        """A hidden client->server channel (a named ``TextInput`` the page writes to)."""
        return self._keep(TextInput(value="", visible=False, name=name))


class StatusLine:
    """A transient message in an HTML bridge that quietly reverts to "Ready".

    Each message replaces the inner node, so the page's CSS animation restarts and
    re-catches the eye; ``level`` (info/warn/danger/success) drives the semantic colour and
    the pop/glow flash. A token guards against an earlier message's revert wiping out a
    newer one. ``undo`` is a removal token (:mod:`..workspace.removal`): the message then
    offers an Undo button carrying it.
    """
    LEVELS = {"info", "warn", "danger", "success"}

    def __init__(self, bridge, doc):
        self.bridge = bridge
        self.doc = doc
        self._token = itertools.count()
        self._last = -1

    def notify(self, message, level="warn", revert_ms=4000, undo=None):
        lvl = level if level in self.LEVELS else "info"
        n = next(self._token)
        self._last = n
        button = ("" if undo is None else
                  f'<button class="status-action" data-undo="{html.escape(str(undo))}">'
                  'Undo</button>')
        # trailing comment keeps the pushed string unique, so an identical repeat
        # message still fires the bridge's change event and replays the flash
        self.bridge.value = (f'<span class="status-msg status-{lvl}">'
                             f'{html.escape(str(message))}{button}</span><!--{n}-->')

        def restore():
            if self._last == n:
                self.bridge.value = "Ready"

        try:
            self.doc.add_timeout_callback(restore, revert_ms)
        except Exception:
            pass


def pulse(bridge, state):
    """Set a state bridge's value to ``"<state>:<nonce>"``.

    The trailing nonce guarantees a change event even when ``state`` is unchanged.
    ``class_bridge``/``check_bridge`` react to the part before the ``:``.
    """
    bridge.value = f"{state}:{next(_nonce)}"


def toggle(bridge, on):
    """Pulse a boolean class/check bridge ``"on"`` or ``"off"``."""
    pulse(bridge, "on" if on else "off")


def parse_action(value):
    """Split an action-channel value ``"<kind> <arg> <nonce>"`` into ``(kind, arg)``.

    The inverse of the page-side ``setv({value: kind + ' ' + arg + ' ' + Date.now()})``
    convention; the trailing nonce (which only forces a change event) is dropped.
    """
    kind, _, rest = value.partition(" ")
    return kind, rest.rsplit(" ", 1)[0]
