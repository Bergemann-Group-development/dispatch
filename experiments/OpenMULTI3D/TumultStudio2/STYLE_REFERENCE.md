# UI Style Reference

A portable guide to the visual language used in this app, so it can be reproduced in a
project with entirely different functionality. It is **app-agnostic** — nothing here
depends on photobooks. It is a compact, dark, dense "pro-tool" aesthetic: think
desktop creative software (editors, DAWs, IDEs) rendered in the browser with plain
HTML/CSS, no framework, no CSS preprocessor.

---

## 1. Design Principles

1. **Dark, low-contrast surfaces.** A near-black background with a small set of slightly
   lighter grays for raised surfaces. Never pure black, never bright white text.
2. **Single accent color.** One blue carries every interactive/active/focus state.
   Color is used sparingly and meaningfully (selection, focus, primary action).
3. **Dense and compact.** Small fonts (10–13px), tight padding (3–10px), thin 1px
   borders. The UI maximizes working area; chrome is quiet.
4. **Flat, with depth only where it floats.** Panels are flat. Shadows appear *only* on
   things that overlay content (modals, canvas-overlay dialogs, the work surface).
5. **Quiet by default, loud on interaction.** Controls sit in muted gray; hover lifts
   them slightly; active/selected states switch to the accent.
6. **Semantic status colors** for success / warning / danger / info, each used both as a
   solid (badges, buttons) and as a translucent tinted background (banners, toasts).

---

## 2. Design Tokens

Defined once as CSS custom properties on `:root`. Copy this block verbatim as the
foundation; everything else references these.

```css
:root {
  /* Layout dimensions (rename/adapt per app) */
  --sidebar-w: 220px;
  --footer-h: 110px;
  --toolbar-h: 40px;

  /* Color palette */
  --bg:          #1e1e1e;   /* app background — deepest layer */
  --surface:     #2a2a2a;   /* raised panels: toolbar, sidebars, modals, footer */
  --border:      #3a3a3a;   /* hairlines, dividers, default control fill */
  --accent:      #4a90e2;   /* the one blue: selection, focus, primary, active */
  --text:        #e0e0e0;   /* primary text (off-white) */
  --text-muted:  #888;      /* secondary text, labels, inactive */
  --canvas-bg:   #3c3c3c;   /* the work surface behind content */
}
```

### Color scale (depth order, dark → light)

| Token / value | Use |
|---|---|
| `#1a1a1a` | recessed wells (preview backings, inset areas) — deeper than `--bg` |
| `--bg` `#1e1e1e` | app background; also the fill for inputs/selects |
| `#242424` | inactive tab background |
| `--surface` `#2a2a2a` | every raised panel and modal |
| `#333` / `--border` `#3a3a3a` | hover fills, hairline borders |
| `#444` | empty thumbnail placeholders |
| `#4a4a4a` / `#505050` | button hover, active toggle fill |
| `#555` / `#606060` / `#666` | control border on hover/emphasis |
| `--text-muted` `#888` | muted text |
| `--text` `#e0e0e0` | primary text |

> Mental model: **every step up in lightness = one step closer to the user.** Background
> `#1e1e1e` → surface `#2a2a2a` → border `#3a3a3a` → hover `#4a4a4a`.

### Accent & semantic colors

| Purpose | Color | Notes |
|---|---|---|
| Accent (blue) | `#4a90e2` | hover-brighter variant `#5a9fe8` / `#5a9fe8` |
| Accent tint (bg) | `rgba(74,144,226,0.1)` – `0.15` | active nav item, selected hover |
| Success (green) | `#22c55e` solid · `#50c878` / `#5fbf63` text · `#2e7d32` border | badges, paid status |
| Warning (amber) | `#e0a030` text · `#e8a020` border · tint `rgba(255,160,0,0.12)` | notes, pending |
| Danger (red) | `#c0392b` solid · `#e05c5c` / `#f06060` text | destructive actions, errors |
| Info (blue toast) | `#80c8f0` text on `#1e2e3a` bg, `#2e5070` border | informational banner |
| Warning banner | `#f0c060` text on `#3a2800` bg, `#7a5000` border | high-contrast attention banner |
| Inline-edit accent | `#34c9a0` (teal) | text-editing caret/outline only |

Pattern for semantic states: **translucent tint background + matching solid border +
matching (brighter) text color.** Example:

```css
.status-pass { background: rgba(80,200,120,0.1); color: #50c878; border: 1px solid #50c878; }
.status-fail { background: rgba(220,60,60,0.12); color: #f06060; border: 1px solid #f06060; }
.status-warn { background: rgba(220,140,0,0.12); color: #f0a030; border: 1px solid #f0a030; }
```

---

## 3. Typography

```css
html, body {
  font: 13px/1.4 system-ui, sans-serif;
  color: var(--text);
  background: var(--bg);
}
```

- **Font family:** `system-ui, sans-serif` — native OS font, no web font for UI.
- **Base size:** `13px`, line-height `1.4`.
- **Monospace** (`monospace`) is reserved for code, logs, and test output.

### Type scale

| Role | Size | Weight | Notes |
|---|---|---|---|
| Big screen title | 30px | normal | start screen hero |
| Hero / overlay title | 18–20px | 600 | loading card, mobile |
| Modal `<h1>`/`<h2>` | 14–15px | 600 | dialog titles |
| App title | 15px | 600 | toolbar brand |
| Card heading | 16px | 600 | |
| Body / inputs | 12–13px | normal | the workhorse size |
| Labels, hints | 10–11px | normal | often muted |
| **Section headers** | **11px** | **600** | `text-transform: uppercase`; color `--text-muted` |
| Tiny captions | 9–10px | normal | thumbnail labels, counts |

### Signature: the uppercase section header

The most recognizable typographic motif. Used for every panel/section title:

```css
.section-header {
  font-size: 11px;
  font-weight: 600;
  text-transform: uppercase;
  color: var(--text-muted);
  letter-spacing: 0.04em;   /* used on some variants */
}
```

---

## 4. Spacing, Radius, Borders

- **Borders:** almost always `1px solid var(--border)`. This single hairline separates
  every region (toolbar bottom, sidebar sides, footer top, modal edges, input outlines).
- **Border radius scale:**
  - `2–3px` — small controls, inputs, chips, thumbnails (the default)
  - `4px` — buttons, standard
  - `6px` — floating overlays / canvas tool palettes
  - `8px` — modals
  - `12px` — large feature cards (start screen, mobile landing)
  - `50%` / `9px` / `10px` — circular badges, pill toggles, status pills
- **Padding rhythm:** controls `3–5px` vertical / `5–10px` horizontal; panels `8–10px`;
  modals `12–24px`. Gaps between items typically `4–8px`.
- **Box-sizing:** global reset `*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }`.

### Shadows (only for floating layers)

```css
--shadow-work:    0 4px 24px rgba(0,0,0,0.5);   /* the content surface */
--shadow-overlay: 0 4px 16px rgba(0,0,0,0.5);   /* canvas-overlay dialogs */
--shadow-modal:   0 8px 32px rgba(0,0,0,0.6);   /* standard modals */
--shadow-modal-lg:0 12px 48px rgba(0,0,0,0.7);  /* large modals */
--shadow-tooltip: 0 2px 8px rgba(0,0,0,0.5);    /* tooltips */
```

Flat panels (sidebars, toolbar, footer) get **no shadow** — they're separated by borders.

---

## 5. Layout Shell

The classic app frame: full-height flex column, fixed-height toolbar + footer, flexible
middle row with sidebars flanking a centered work surface.

```css
html, body { height: 100%; overflow: hidden; }

#app      { display: flex; flex-direction: column; height: 100vh; }

#toolbar  { height: var(--toolbar-h); background: var(--surface);
            border-bottom: 1px solid var(--border);
            display: flex; align-items: center; padding: 0 12px; gap: 12px;
            flex-shrink: 0; }

#workspace{ flex: 1; display: flex; overflow: hidden; }   /* middle row */

aside     { width: var(--sidebar-w); background: var(--surface);
            border-right: 1px solid var(--border);
            display: flex; flex-direction: column; overflow: hidden; flex-shrink: 0; }
/* right sidebar mirrors the border to the left edge */
.sidebar-right { border-right: none; border-left: 1px solid var(--border); }

#main     { flex: 1; overflow: hidden; display: flex;
            align-items: center; justify-content: center; background: var(--canvas-bg); }

#footer   { height: var(--footer-h); background: var(--surface);
            border-top: 1px solid var(--border);
            display: flex; align-items: center; padding: 0 8px; gap: 6px; flex-shrink: 0; }
```

Key idioms:
- `flex-shrink: 0` on toolbar/footer/sidebars so only the center flexes.
- `overflow: hidden` on shell containers; scrolling lives on inner content areas
  (`overflow-y: auto` with `flex: 1; min-height: 0`).
- The work surface (`--canvas-bg`, lighter than `--bg`) is visually distinct from chrome.
- The "document" sitting on the surface gets `box-shadow: 0 4px 24px rgba(0,0,0,0.5)` to
  float above it.

### Toolbar internals
- `.app-title` left (15px/600), action groups with `gap: 6px`, and a
  `.toolbar-right { margin-left: auto }` pushing controls to the far right.

---

## 6. Components

### 6.1 Buttons

Base button = muted gray, hover lightens, `.active`/`.primary` = accent.

```css
button {
  background: var(--border);
  color: var(--text);
  border: 1px solid #555;
  border-radius: 4px;
  padding: 4px 10px;
  font-size: 12px;
  cursor: pointer;
  transition: background 0.1s;
}
button:hover         { background: #4a4a4a; }
button.active        { background: var(--accent); border-color: var(--accent); }
button.primary       { background: var(--accent); border-color: var(--accent); color: #fff; }
button.primary:hover { background: #5a9fe8; }
```

**Secondary / ghost button** (very common — used in dialogs, cards): transparent or
surface fill, border lightens to accent on hover, text turns accent:

```css
.btn-ghost {
  background: var(--surface);          /* or transparent */
  border: 1px solid var(--border);
  border-radius: 4px;
  color: var(--text);
  padding: 5px 14px;
  cursor: pointer;
}
.btn-ghost:hover { border-color: var(--accent); color: var(--accent); }
```

> **The signature hover.** Across cards, nav items, confirm buttons, and icon buttons the
> recurring hover effect is *not* a background change but **border + text both turning to
> the accent color.** This is the app's most consistent interaction tell.

**Danger / confirm variants** (filled, hover dims via opacity):

```css
.btn-danger  { background: #c0392b;      border-color: #c0392b;      color: #fff; }
.btn-confirm { background: var(--accent); border-color: var(--accent); color: #fff; }
.btn-danger:hover, .btn-confirm:hover { opacity: 0.85; }
```

**Small icon buttons** — square, centered glyph, ghost-hover:

```css
.icon-btn {
  height: 22px; width: 22px; padding: 0;
  display: flex; align-items: center; justify-content: center;
  background: var(--bg);
  border: 1px solid var(--border);
  border-radius: 3px;
  color: var(--text-muted);
  cursor: pointer;
}
.icon-btn:hover { border-color: var(--accent); color: var(--accent); }
.icon-btn:disabled { opacity: 0.3; cursor: default; }
```

### 6.2 Form fields (input / select / textarea)

All text/number/select/textarea inputs share one recipe: dark `--bg` fill, hairline
border, accent border on focus, no focus ring.

```css
input[type="number"], input[type="text"], select, textarea {
  width: 100%;
  background: var(--bg);
  border: 1px solid var(--border);
  border-radius: 3px;
  color: var(--text);
  padding: 3px 5px;
  font-size: 12px;
}
input:focus, select:focus, textarea:focus {
  outline: none;
  border-color: var(--accent);
}
```

- **Field group:** label stacked above input, `display: flex; flex-direction: column; gap: 2px`.
  Labels are `10px`, `--text-muted`.
- **Color inputs:** `width: 100%` (or fixed e.g. 36–50px), `height: ~22–28px`, small padding.
- **Two-column field grid:** `display: grid; grid-template-columns: 1fr 1fr; gap: 4px;`
  with a `.full-width { grid-column: 1 / -1; }` escape hatch.
- **"Mixed value" state** (multi-select disagreement): dashed dimmed border on color
  swatches, italic muted placeholder on number inputs — a nice pattern for any
  multi-target editor.

### 6.3 Toggle switch (pill)

```css
.switch { position: relative; display: inline-block; width: 32px; height: 18px; cursor: pointer; }
.switch input { position: absolute; opacity: 0; width: 0; height: 0; }
.switch-track { position: absolute; inset: 0; background: var(--border);
                border-radius: 9px; transition: background 0.15s ease; }
.switch-track::before {
  content: ""; position: absolute; width: 14px; height: 14px; left: 2px; top: 2px;
  background: #fff; border-radius: 50%; transition: transform 0.15s ease;
}
.switch input:checked + .switch-track          { background: var(--accent); }
.switch input:checked + .switch-track::before   { transform: translateX(14px); }
```

### 6.4 Segmented / toggle-button groups

Rows of equal-width buttons; selected one fills with accent (or a neutral `#505050` for
secondary segmentation), `gap: 2px`.

```css
.seg-row { display: flex; gap: 2px; height: 26px; }
.seg-btn {
  flex: 1; padding: 0;
  background: var(--bg);
  border: 1px solid var(--border);
  border-radius: 3px;
  color: var(--text-muted);
  cursor: pointer;
}
.seg-btn.active        { background: var(--accent); border-color: var(--accent); color: #fff; }
.seg-btn:hover:not(.active) { border-color: #666; color: var(--text); }
/* neutral variant for non-primary segmentation: */
.seg-btn--neutral.active { background: #505050; color: var(--text); border-color: #606060; }
```

### 6.5 Tabs

Two flavors:

**Classic tab strip** (active tab merges into the panel below by matching its bg and
hiding its bottom border):

```css
.tab {
  flex: 1; padding: 8px 10px;
  font-size: 11px; font-weight: 600; text-transform: uppercase;
  color: var(--text-muted); background: #242424;
  border: none; border-right: 1px solid var(--border); border-bottom: 1px solid var(--border);
  cursor: pointer;
}
.tab.active { color: var(--text); background: var(--surface); border-bottom-color: var(--surface); }
```

**Underline tabs** (for dialogs): borderless, 2px accent underline on the active tab:

```css
.utab { padding: 7px 14px; border: none; border-bottom: 2px solid transparent;
        background: none; color: var(--text-muted); font-size: 13px; cursor: pointer; }
.utab.active { color: var(--text); border-bottom-color: var(--accent); }
```

### 6.6 Selectable items / thumbnails

Selection is shown with a **2px transparent border that turns accent when selected**
(border-box reserves the space so nothing shifts):

```css
.item {
  border: 2px solid transparent;
  border-radius: 3px;
  cursor: pointer;
  transition: border-color 0.1s;
}
.item:hover    { border-color: #666; }
.item.selected { border-color: var(--accent); }
```

A small **corner badge** pattern (e.g. status tick): absolutely positioned circle in the
corner, solid semantic color, white glyph.

```css
.badge {
  position: absolute; top: 3px; right: 3px;
  width: 16px; height: 16px; border-radius: 50%;
  background: #22c55e; color: #fff;
  font-size: 10px; font-weight: 700; line-height: 16px; text-align: center;
}
```

### 6.7 Status pills

```css
.pill {
  font-size: 11px; padding: 2px 8px; border-radius: 10px;
  background: var(--surface); border: 1px solid var(--border);
  white-space: nowrap;
}
.pill-ok      { border-color: #2e7d32; color: #5fbf63; }
.pill-pending { border-color: #e0a030; color: #e0a030; }
.pill-fail    { border-color: #c0392b; color: #e05c5c; }
```

### 6.8 Progress bars

A shared track/fill pattern — track is `--border`, fill is `--accent`:

```css
.progress-track { background: var(--border); overflow: hidden; border-radius: 3px; height: 6px; }
.progress-bar   { width: 0%; height: 100%; background: var(--accent);
                  border-radius: 3px; transition: width 0.15s ease-out; }
```

Indeterminate variant: a 25%-wide bar sliding via keyframes (`translateX(-100%) → 400%`).

---

## 7. Overlays

### 7.1 Modals / dialogs

Built on the native `<dialog>` element. Shared chrome:

```css
.modal {
  background: var(--surface);
  color: var(--text);
  border: 1px solid var(--border);
  border-radius: 8px;
  box-shadow: 0 8px 32px rgba(0,0,0,0.6);
  position: fixed; top: 50%; left: 50%; translate: -50% -50%;
}
.modal::backdrop { background: rgba(0,0,0,0.5); }
```

- **Sizing:** responsive caps via `width: min(480px, 92vw)` and
  `max-height: min(560px, 88vh)`. Pick a px ceiling per dialog; the `vw`/`vh` keeps it
  safe on small screens.
- **Structure:** header (with bottom border) / scrollable content (`overflow-y: auto;
  flex: 1; min-height: 0`) / action footer (with top border). Header & footer
  `flex-shrink: 0`. Set `padding: 0` on the dialog and pad each region instead.
- **Header:** `display: flex; justify-content: space-between; align-items: center;
  padding: 12px 16px 10px; border-bottom: 1px solid var(--border);` title 14px/600.
- **Action row:** `display: flex; justify-content: flex-end; gap: 8px;` — secondary
  (ghost) button(s) then primary/danger on the right.
- **Sizes seen:** confirm `min(300px,92vw)`, standard `min(480px,92vw)`, settings
  `min-width:300px`, two-column settings `min(720px,92vw)`, docs `min(960px,92vw)` ×
  `min(680px,88vh)`.

### 7.2 Canvas-overlay (floating) dialogs & tool palettes

Small panels that float over the work surface — **translucent dark + backdrop blur**,
the distinguishing treatment vs. opaque modals:

```css
.tool-dialog, .floating-palette {
  position: absolute;
  background: rgba(30,30,30,0.85);     /* dialogs use 0.92 */
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 4px 6px;                    /* dialogs ~8–10px */
  backdrop-filter: blur(4px);
  z-index: 10;
}
```

Used for zoom controls, tool palettes, and contextual mini-dialogs anchored to corners
of the work area (`top: 22px; left/right: …`).

### 7.3 Tooltips

```css
.tooltip {
  position: fixed; z-index: 9999;
  background: #2a2a2a; color: var(--text);
  border: 1px solid var(--border); border-radius: 4px;
  padding: 5px 9px; font-size: 12px;
  pointer-events: none; white-space: nowrap;
  box-shadow: 0 2px 8px rgba(0,0,0,0.5);
}
```

### 7.4 Banners (inline, full-width)

Slim attention strips with `flex; align-items: center; gap: 8px; padding: 6px 10px;
font-size: 12px;` and a `margin-left:auto` dismiss `×` button. Colored by intent using
the tint+border+text recipe — e.g. warning `bg #3a2800 / border #7a5000 / text #f0c060`,
info `bg #1e2e3a / border #2e5070 / text #80c8f0`.

### 7.5 Full-screen overlays (loading, preview)

`position: fixed; inset: 0;` with `--bg` (or near-black `#0d0d0d` for an immersive
preview) and centered flex content. Use a high `z-index` (9999 for loaders).

---

## 8. Content typography (rich text / docs)

For long-form rendered content (markdown, help text) the prose styles are:

```css
.prose            { line-height: 1.65; }
.prose h1         { font-size: 20px; font-weight: 600; margin-bottom: 16px; }
.prose h2         { font-size: 15px; font-weight: 600; margin: 24px 0 10px;
                    border-bottom: 1px solid var(--border); padding-bottom: 4px; }
.prose h3         { font-size: 13px; font-weight: 600; margin: 18px 0 8px; }
.prose p          { margin-bottom: 12px; }
.prose a          { color: var(--accent); text-decoration: none; }
.prose a:hover    { text-decoration: underline; }
.prose code       { background: var(--bg); border: 1px solid var(--border);
                    border-radius: 3px; padding: 1px 5px; font: 11px monospace; }
.prose pre        { background: var(--bg); border: 1px solid var(--border);
                    border-radius: 4px; padding: 12px; overflow-x: auto; }
.prose blockquote { border-left: 3px solid var(--accent); padding: 6px 14px;
                    background: rgba(74,144,226,0.07); color: var(--text-muted); }
.prose th         { background: var(--bg); border: 1px solid var(--border);
                    padding: 6px 10px; font-weight: 600; text-align: left; }
.prose td         { border: 1px solid var(--border); padding: 5px 10px; }
```

---

## 9. Iconography

- **Library:** [Tabler Icons](https://tabler.io/icons) via web-font CDN:
  ```html
  <link rel="stylesheet"
        href="https://cdn.jsdelivr.net/npm/@tabler/icons-webfont@latest/dist/tabler-icons.min.css" />
  ```
  Usage: `<i class="ti ti-folder-open" aria-hidden="true"></i>`.
- Icons inherit `currentColor`, so they pick up the muted/accent text color of their
  control automatically.
- A spin utility for loading glyphs:
  ```css
  @keyframes ti-spin { to { transform: rotate(360deg); } }
  .ti-spin { animation: ti-spin 1s linear infinite; }
  ```

---

## 10. Motion

Animations are minimal and fast — they confirm interaction without drawing attention.

- **Hover/state transitions:** `transition: background 0.1s;` or `border-color 0.1s;`
  (~100ms).
- **Toggle switches:** `0.15s ease`.
- **Progress fills:** `width 0.15s ease-out`.
- No easing flourishes, no large movement, no entrance animations on panels.

---

## 11. Responsive / Mobile

This is a desktop-first tool. Below a breakpoint it hides the app entirely and shows a
single centered "use a larger screen" card:

```css
@media (max-width: 1023px) {
  #app, #start-screen { display: none !important; }
  #mobile-landing { display: flex; align-items: center; justify-content: center;
                    min-height: 100vh; padding: 24px; }
}
```

Feature cards (start screen / mobile) use the larger `12px` radius, `--surface` fill,
centered column layout, accent-colored hero icon (30–48px), and the ghost-button hover.

---

## 12. Quick-start checklist for a new app

1. Drop in the `:root` token block (§2) and the global box-sizing reset.
2. Set `body` font to `13px/1.4 system-ui` with `--bg` / `--text`.
3. Build the shell: flex-column `#app`, fixed toolbar/footer, flex `#workspace` with
   `aside` sidebars around a `--canvas-bg` `#main` (§5).
4. Add the base `button`, form-field, and `.section-header` recipes (§3, §6.1, §6.2).
5. Reach for the semantic tint+border+text trio for any status UI (§2).
6. Use native `<dialog>` + the `.modal` chrome for dialogs; translucent+blur for things
   floating over the work surface (§7).
7. Pull icons from Tabler web-font; let them inherit `currentColor`.
8. Keep transitions ~100–150ms; shadows only on floating layers.

**The five things that make it look like *this* app:** (1) the gray depth ladder
`#1e1e1e→#2a2a2a→#3a3a3a`, (2) the single blue accent `#4a90e2`, (3) uppercase 11px/600
muted section headers, (4) the border+text→accent hover, and (5) translucent
blurred panels floating over a lighter work surface.
