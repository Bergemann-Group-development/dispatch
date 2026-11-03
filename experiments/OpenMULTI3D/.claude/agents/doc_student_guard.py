#!/usr/bin/env python3
"""PreToolUse guard for the doc-student agent.

The agent works in a fresh clone whose root holds a `.student-sandbox` marker. It may
read the user documentation, the inputs and outputs the documentation points to and the
files it created itself; it may not read the source code or leave the sandbox. Shell use
is limited to the commands the installation guide and quickstart type at a prompt.

Exit code 2 blocks the call; stderr is returned to the agent as the reason.
"""

import json
import re
import shlex
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

MARKER = ".student-sandbox"
PROJECT = Path("experiments/OpenMULTI3D")

# User documentation, relative to the project directory (prefixes end in "/").
DOCS = (
    "README_INSTALL.md",
    "QUICKSTART.md",
    "docs/FAQ.md",
    "docs/GLOSSARY.md",
    "docs/GUIWorkflow.md",
    "docs/manual/",
    "HOP/README.md",
)
# Tracked inputs the documentation tells the reader to open.
DATA = ("input/", "data/")
STUDENT_SUFFIXES = {".py", ".ipynb", ".txt", ".md", ".nml", ".json", ".log", ".csv",
                    ".dat", ".png", ".jpg", ".pdf", ".svg"}
WRITE_SUFFIXES = {".py", ".ipynb", ".txt", ".md", ".nml", ".json", ".csv", ".dat"}

FIRST_WORDS = {
    "cd", "ls", "pwd", "echo", "export", "which", "source", ".", "deactivate",
    "python", "python3", "pip", "pip3", "make", "./dispatch.x", "git", "gfortran",
    "perl", "brew", "jupyter", "bokeh", "mkdir", "cp", "mv", "rm", "kill", "pkill",
    "curl", "conda", "cat", "head", "tail", "less", "more", "grep", "wc", "sort",
    "file", "du", "df", "uname", "sw_vers", "sleep", "true",
}
FILE_READERS = {"cat", "head", "tail", "less", "more", "grep", "wc", "sort", "file"}
GIT_SUBCOMMANDS = {"clone", "status", "log", "pull", "fetch", "rev-parse", "branch",
                   "checkout", "switch", "remote"}
SHELL_SYNTAX = re.compile(r"<<|\b(for|while|until|function|case)\b.*\b(do|in)\b|\$\(|`")
SOURCE_PEEK = re.compile(r"getsource|inspect\.getfile|\.f90|\.i90|__code__|dis\.dis")


def deny(reason):
    print(f"Blocked for the student simulation: {reason}", file=sys.stderr)
    sys.exit(2)


def sandbox_of(path):
    for parent in (path, *path.parents):
        if (parent / MARKER).is_file():
            return parent
    return None


@lru_cache(maxsize=None)
def tracked(root):
    out = subprocess.run(["git", "-C", str(root), "ls-files"], capture_output=True,
                         text=True, check=True).stdout
    return frozenset(out.splitlines())


def locate(raw, cwd):
    path = Path(raw).expanduser()
    path = (path if path.is_absolute() else cwd / path).resolve()
    root = sandbox_of(path)
    if root is None:
        deny(f"{raw} lies outside the sandbox clone.")
    return path, root


def relative(path, root):
    """Path relative to the project directory, or None outside it."""
    try:
        return path.relative_to(root / PROJECT).as_posix()
    except ValueError:
        return None


def is_doc(rel):
    return any(rel == d or (d.endswith("/") and rel.startswith(d)) for d in DOCS)


def student_file(path, root, rel, suffixes):
    in_git = path.relative_to(root).as_posix() in tracked(root)
    return (not in_git and rel is not None and not rel.startswith(".venv/")
            and path.suffix in suffixes)


def check_read(raw, cwd):
    path, root = locate(raw, cwd)
    if path == root / "README.md":
        return
    rel = relative(path, root)
    if rel is None:
        deny(f"{raw} is not part of OpenMULTI3D's user documentation.")
    if path.is_dir():
        if rel == "docs" or rel.startswith(("docs/",) + DATA):
            return
        deny(f"{rel or '.'}: only the docs/, input/ and data/ folders may be searched.")
    if is_doc(rel) or rel.startswith(DATA) or (
            "/" not in rel and path.suffix == ".nml"):
        return
    if student_file(path, root, rel, STUDENT_SUFFIXES):
        return
    deny(f"{rel} is source code or developer material, not user documentation. "
         "Record what you wanted to learn from it as an obstacle.")


def check_write(raw, cwd):
    path, root = locate(raw, cwd)
    rel = relative(path, root)
    if rel is None:
        deny(f"write {raw} inside the OpenMULTI3D project directory.")
    if path.suffix == ".nml" and "/" not in rel:
        return
    if student_file(path, root, rel, WRITE_SUFFIXES):
        return
    deny(f"{rel}: write only your own Python scripts, notebooks, namelists and notes, "
         "never shell scripts or tracked source files.")


def check_python(text):
    if SOURCE_PEEK.search(text):
        deny("this Python code inspects the package source. Use help() or the "
             "documentation instead.")


def segments(command):
    """The simple commands of a command line, split at operators outside quotes."""
    lexer = shlex.shlex(command, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    words = []
    try:
        for token in lexer:
            if token and set(token) <= set("&|;()"):
                yield words
                words = []
            elif not (token and set(token) <= set("<>")):
                words.append(token)
    except ValueError:
        deny("the command could not be parsed.")
    yield words


def check_bash(command, cwd):
    if SHELL_SYNTAX.search(command):
        deny("shell scripting (loops, heredocs, command substitution) is outside the "
             "student's skills. Type the documented commands, or write Python.")
    check_python(command)
    for words in segments(command):
        while words and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", words[0]):
            words.pop(0)
        if not words:
            continue
        head, args = words[0], words[1:]
        name = head if head.startswith("./") else Path(head).name
        if name not in FIRST_WORDS and not head.endswith("/bin/python"):
            deny(f"'{head}' is not a command the documentation uses.")
        if name in ("bash", "sh", "zsh") or any(a.endswith(".sh") for a in words):
            deny("shell scripts are not allowed.")
        if name == "brew" and args[:1] not in (["list"], ["--version"], ["--prefix"]):
            deny("do not install system packages here. Check whether the tool is "
                 "present (e.g. `gfortran --version`) and record the step.")
        if name == "git" and args[:1] and args[0] not in GIT_SUBCOMMANDS:
            deny(f"'git {args[0]}' is not needed to follow the documentation.")
        if name == "curl" and not any("localhost" in a or "127.0.0.1" in a for a in args):
            deny("curl is only for checking that the local GUI server answers.")
        if name == "./dispatch.x" and any("3d" in a.lower() for a in args):
            deny("3D runs are not allowed on this machine. Record the 3D example as "
                 "not tested.")
        if name in ("python", "python3") or head.endswith("/bin/python"):
            scripts = [a for a in args if a.endswith(".py")]
            if args[:1] not in (["-c"], ["-m"]) and scripts:
                check_read(scripts[0], cwd)
                check_python(Path(cwd, scripts[0]).read_text(errors="ignore")
                             if Path(cwd, scripts[0]).is_file() else "")
        if name in FILE_READERS | {"cp", "mv"}:
            operands = [a for a in args if not a.startswith("-")]
            for a in operands[1:] if name == "grep" else operands:
                if Path(cwd, a).expanduser().is_file() or "/" in a or "." in a:
                    check_read(a, cwd)
        if name == "rm":
            for a in args:
                if not a.startswith("-"):
                    path, root = locate(a, cwd)
                    rel = relative(path, root) or ""
                    if not (rel.startswith("data/") or student_file(
                            path, root, rel, STUDENT_SUFFIXES) or "*" in a):
                        deny(f"rm {a}: remove only your own files and run outputs.")


def main():
    event = json.load(sys.stdin)
    tool, args = event.get("tool_name", ""), event.get("tool_input", {})
    cwd = Path(event.get("cwd", ".")).resolve()
    if sandbox_of(cwd) is None:
        deny(f"the working directory {cwd} is outside the sandbox clone. cd into it.")
    if tool == "Read":
        check_read(args["file_path"], cwd)
    elif tool in ("Write", "Edit", "NotebookEdit"):
        target = args.get("file_path") or args.get("notebook_path")
        check_write(target, cwd)
        check_python(args.get("content", "") + args.get("new_string", "")
                     + args.get("new_source", ""))
    elif tool == "Grep":
        raw = args.get("path") or deny("give Grep a documentation file or folder.")
        check_read(raw, cwd)
    elif tool == "Glob":
        locate(args.get("path") or ".", cwd)
    elif tool == "Bash":
        check_bash(args.get("command", ""), cwd)
    sys.exit(0)


if __name__ == "__main__":
    main()
