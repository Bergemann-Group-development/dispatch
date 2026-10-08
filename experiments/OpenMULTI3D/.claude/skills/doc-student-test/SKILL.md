---
name: doc-student-test
description: Test whether the OpenMULTI3D user documentation lets a new student install and run the code, by launching the doc-student agent in a fresh sandbox clone and reviewing its journal. Use when asked to test the documentation's accessibility, the installation guide or the quickstart.
---

# Documentation test with the doc-student agent

The `doc-student` agent (`.claude/agents/doc-student.md`) follows the user documentation
as a student with basic Python would, without reading source code. Its guard hook
(`.claude/agents/doc_student_guard.py`) only allows actions inside a directory tree whose
root holds a `.student-sandbox` marker. The agent changes no tracked code, so it runs in a
plain clone and needs none of the worktree and branch procedure of `AGENTS.md`.

## 1. Prepare the sandbox

The sandbox is a clone of a commit, the way a student receives the repository. Commit
documentation edits that are to be tested first; uncommitted changes are not in it.

```sh
base=$(git rev-parse HEAD)
sandbox=<scratchpad>/doc-student-$(date +%Y%m%d-%H%M)
git clone -q "$(git rev-parse --show-toplevel)" "$sandbox"
git -C "$sandbox" checkout -q --detach "$base"
touch "$sandbox/.student-sandbox"
echo .student-sandbox >> "$sandbox/.git/info/exclude"
```

Use the session scratchpad, or a path the user names. The clone holds the tracked input
data and needs a few hundred MB; the build and `.venv` add about as much again.

## 2. Launch

Launch `doc-student` in the background with an assignment that names:

- the sandbox root and the base commit;
- the scope, by default `README_INSTALL.md` then `QUICKSTART.md` including its optional
  sections except the 3D snapshot; or a narrower one such as a single FAQ answer;
- anything the user wants varied, e.g. "you have never used a virtual environment".

Expect 30 to 60 minutes; the native build alone takes several minutes. `model: sonnet` in
the agent file is the default; `haiku` gets stuck more readily and gives a lower bound on
accessibility, at the cost of more obstacles that are its own mistakes.

## 3. Review

Read `<sandbox>/experiments/OpenMULTI3D/student_journal.md` and the agent's report. For
each obstacle decide its cause before reporting it:

- **documentation**: wrong, missing, out of order or unexplained;
- **code**: the documentation is right and the program misbehaves;
- **environment**: specific to this machine (conda on `PATH`, installed compilers,
  `OMP_NUM_THREADS=1` in the agent's shell);
- **agent**: the simulated student misread a clear instruction, or a harness effect.

Reproduce any blocker or deviation that is not plainly a documentation issue, in the
sandbox, before reporting it. Report to the user grouped by cause, with file and section,
and list the steps that were not tested. Remove the sandbox once the user has the report.
