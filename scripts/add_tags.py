#!/usr/bin/env python3
# --------------------------------------------------------------------------
# add_tags.py
# Oral microbiome age clock: cross-study validation of salivary 16S models
#
# Francisco M. Codoner
# Instituto de Agroquimica y Tecnologia de Alimentos (IATA-CSIC), Valencia, Spain
# ORCID 0000-0003-3671-8148
#
# Repository: https://github.com/fcodoner/oral-age-clock
# Released under the MIT License; see LICENSE.
# --------------------------------------------------------------------------
"""
add_tags.py - add metadata labels to one of the *_TAGS lists in
harvest_oral_age_metadata.py, safely, without hand-editing.

  python scripts/add_tags.py SITE_TAGS "environment (material)" "tissue"
  python scripts/add_tags.py SEX_TAGS "host_sex" "gender"

Makes a .bak backup, keeps the existing entries, skips duplicates,
rewrites the list on one properly indented line, and checks the file compiles.
"""
import re, shutil, sys, py_compile

SCRIPT = "scripts/harvest_oral_age_metadata.py"


def main():
    if len(sys.argv) < 3:
        sys.exit('Usage: python scripts/add_tags.py SITE_TAGS "label one" "label two"')
    name, new_tags = sys.argv[1], sys.argv[2:]
    try:
        src = open(SCRIPT).read()
    except FileNotFoundError:
        sys.exit(f"Cannot find {SCRIPT}. Are you in the ~/oral-age-clock folder?")

    m = re.search(rf"^{name}\s*=\s*\[(.*?)\]", src, re.S | re.M)
    if not m:
        lists = re.findall(r"^([A-Z_]+_TAGS)\s*=", src, re.M)
        sys.exit(f"No list called {name}. Available: {', '.join(sorted(set(lists)))}")

    existing = [t.strip().strip('"\'') for t in m.group(1).split(",") if t.strip()]
    added = [t for t in new_tags if t.lower() not in [e.lower() for e in existing]]
    if not added:
        print("Nothing to add - all of those labels are already in the list.")
        return
    merged = existing + added
    body = ", ".join(f'"{t}"' for t in merged)
    new_line = f"{name} = [{body}]"

    shutil.copy(SCRIPT, SCRIPT + ".bak")
    open(SCRIPT, "w").write(src[:m.start()] + new_line + src[m.end():])
    try:
        py_compile.compile(SCRIPT, doraise=True)
    except py_compile.PyCompileError as e:
        sys.exit(f"Edit broke the file:\n{e}\nRestore with:  mv {SCRIPT}.bak {SCRIPT}")
    print(f"Added to {name}: " + ", ".join(added))
    print(f"{name} now has {len(merged)} labels. Backup at {SCRIPT}.bak")


if __name__ == "__main__":
    main()
