"""A tiny server->client bridge: a hidden widget whose value is pushed as
innerHTML into a light-DOM element (so it can be styled with ordinary CSS and
support event delegation, unlike Bokeh's shadow-DOM widgets)."""
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

    def check(self, input_id):
        return self._keep(check_bridge(input_id))

    def action(self, name):
        """A hidden client->server channel (a named ``TextInput`` the page writes to)."""
        return self._keep(TextInput(value="", visible=False, name=name))


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
