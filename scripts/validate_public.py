#!/usr/bin/env python3
"""Fail when material that must remain private appears in the public tree."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FORBIDDEN_NAME_MARKERS = (
    "instructor-only",
    "private",
    "solution",
    "unreleased",
)
FORBIDDEN_NOTEBOOK_TAGS = {
    "instructor-only",
    "private",
    "solution",
}
FORBIDDEN_SENTINELS = (
    "INSTRUCTOR_ONLY",
    "PRIVATE_ONLY",
)
FORBIDDEN_SUFFIXES = {
    ".aux",
    ".bib",
    ".fls",
    ".log",
    ".synctex",
    ".tex",
}
IGNORED_PARTS = {".git", ".ipynb_checkpoints", "__pycache__", "dist"}


def iter_files(root: Path):
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if any(part in IGNORED_PARTS for part in relative.parts):
            continue
        if path.is_symlink():
            yield path, "symbolic links are not allowed in public material"
        elif path.is_file():
            yield path, None


def validate_notebook(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        notebook = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return [f"invalid notebook JSON: {exc}"]

    if notebook.get("nbformat") != 4 or not isinstance(notebook.get("cells"), list):
        errors.append("expected notebook format 4 with a cell list")

    for index, cell in enumerate(notebook.get("cells", [])):
        tags = set(cell.get("metadata", {}).get("tags", []))
        leaked = tags & FORBIDDEN_NOTEBOOK_TAGS
        if leaked:
            errors.append(f"cell {index} retains private tags: {sorted(leaked)}")

    serialized = json.dumps(notebook, ensure_ascii=False)
    for sentinel in FORBIDDEN_SENTINELS:
        if sentinel in serialized:
            errors.append(f"contains private sentinel {sentinel!r}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", default=".")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    errors: list[str] = []

    for path, path_error in iter_files(root):
        relative = path.relative_to(root)
        if path_error:
            errors.append(f"{relative}: {path_error}")
            continue

        lowered_parts = [part.casefold() for part in relative.parts]
        for part in lowered_parts:
            if any(marker in part for marker in FORBIDDEN_NAME_MARKERS):
                errors.append(f"{relative}: forbidden private-looking filename")
                break

        suffixes = {suffix.casefold() for suffix in path.suffixes}
        if suffixes & FORBIDDEN_SUFFIXES:
            errors.append(f"{relative}: authoring source or build by-product")

        if path.suffix.casefold() == ".ipynb":
            errors.extend(f"{relative}: {error}" for error in validate_notebook(path))

    if errors:
        print("Public-content validation failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("Public-content validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

