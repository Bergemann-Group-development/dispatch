---
name: manual-diagrams
description: Draw flowcharts and task diagrams for the OpenMULTI3D manual (docs/manual) as TikZ figures in the house style. Use when a manual section needs a diagram, or when an existing figure (including the old draw.io PDFs in docs/manual/figs) is to be redrawn or restyled.
---

# Diagrams for the manual

Diagrams are TikZ files in `docs/manual/figs/<name>.tex`, pulled into `main.tex` with
`\input`. `figs/passes.tex` (Fig. "Wavelength grids of the two passes") and
`figs/tasks.tex` (Fig. "Tasks of a run") are the reference flowcharts, and
`figs/dispatch_step.tex` (Fig. "Calls of one dispatcher step") the reference call tree;
start from the one closer to the new diagram.

## Content first

Read the manual section the figure belongs to and the code it describes before drawing.
A diagram states facts in compressed form, so each box and arrow needs a source: a
dependency in a `nbor_relations`, a branch in the setup code, a key in a namelist table.
Where the manual and the code disagree, say so to the user instead of drawing either.

## File skeleton

```latex
% <what the figure shows> (Sec. "<section title>").
% Needs tikz with the libraries positioning, arrows.meta, calc and fit.
\begingroup
\newcommand{\kw}[1]{\texttt{#1}}
\begin{adjustbox}{max width=\linewidth, max height=0.62\textheight}
\begin{tikzpicture}[
    font=\footnotesize,
    >={Stealth[length=2mm]},
    box/.style={draw, rounded corners=4pt, align=center, inner sep=4pt, text width=5.5cm},
    arr/.style={->, thick},
    lab/.style={font=\scriptsize, text=black!70, align=left},
    head/.style={font=\small\bfseries, align=center},
    node distance=5mm
  ]
  ...
\end{tikzpicture}
\end{adjustbox}
\endgroup
```

and in `main.tex`, next to the text it illustrates:

```latex
\begin{figure}[!ht]
    \centering
    \input{figs/<name>}
    \caption{<one or two sentences: what is shown and the one point to take away>}\label{fig:<name>}
\end{figure}
```

Write namelist keys with `\kw{...}` (escape `_` and `&`). `\fort` is verbatim and breaks
inside TikZ nodes. Do not name a style `step`: it is a TikZ key.

## Colours

The palette lives in the `main.tex` preamble as the styles `fig grey`, `fig blue`,
`fig orange`, `fig green`, `fig red` and `fig purple`: a pastel fill with a darker border of the same
hue, after the draw.io defaults. Define per-figure styles on top of them
(`grid/.style={box, fig blue}`) and keep the meaning of each hue:

| style | used for |
|---|---|
| `fig grey` | processing steps (window splits and the like); DISPATCH routines in call diagrams |
| `fig blue` | wavelength grids |
| `fig orange` | input data, composition |
| `fig green` | NLTE: iterations, populations, departure coefficients |
| `fig red` | radiative transfer and its results |
| `fig purple` | the OpenMULTI3D task layer (`m3d_patch_t`) in call diagrams |

A new kind of box gets a new palette entry in the preamble, not an ad-hoc colour in the
figure. Plain arrows are black; arrows inside a coloured group may take its line colour
(`figRedLine`).

## Text in a box

The first line or a bold title carries the main statement in black. Everything secondary
(descriptions, parameter lists, file names) goes in `\fignote{...}`, which renders in a
darker shade of the box's hue. A `\\` inside `\fignote` works only in a node with a
`text width`; in a node sized to its content, put each line in its own `\fignote`.

```latex
\node[grid] (bb) {\textbf{bound-bound grid}\\
  \fignote{profile points over each active line (\kw{qnorm}, \kw{nq\_fac})}};
```

## Branches

No diamonds. A branch is a labelled fork: the condition sits centred *on* the stem, the
stem continues to a fork point, and each leg carries its value beside the downward part.

```latex
\coordinate (dc) at ($(atoms.south)+(0,-10mm)$);
\draw[thick] (atoms) -- (dc) node[midway, lab, fill=white, inner sep=1.5pt] {\kw{decouple\_continuum}};
\draw[arr] (dc) -| (bb.north)  node[pos=0.75, right=1pt, lab] {T (default)};
\draw[arr] (dc) -| (all.north) node[pos=0.75, right=1pt, lab] {F};
```

## Arrows

- Every arrow ends on its target box. A box that receives arrows from several branches
  is centred under the fork, not under one branch, and is wide enough for all of them.
- Arrows do not cross where a straight route exists; prefer `(a.south) -- (a.south |- b.north)`
  for vertical drops into a wide box.
- Labels never sit on top of a line or a box edge. Give arrows enough length for their
  label, or put the label on a white fill interrupting the line.

## Check the rendering

Build the manual and look at the page before reporting. From `docs/manual`:

```sh
latexmk -pdf -shell-escape -interaction=nonstopmode main.tex | grep -E "^!|Output written"
python3 - <<'EOF'
import fitz
for i, p in enumerate(fitz.open('main.pdf')):
    if '<words from the caption>' in p.get_text():
        print(i + 1); p.get_pixmap(dpi=110).save('<scratchpad>/page.png')
EOF
```

Read the PNG and check: arrows connect, nothing crosses needlessly, no label overlaps a
line or edge, the float stays near its section (a figure that drifts to the end is too
tall), and the text is legible at the scaled size. Iterate until it is clean.

The parent repository's `.gitignore` excludes much of `docs/`; check `git status` shows a
new `figs/*.tex` before committing.
