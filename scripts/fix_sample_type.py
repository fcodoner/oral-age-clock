#!/usr/bin/env python3
# --------------------------------------------------------------------------
# fix_sample_type.py
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
fix_sample_type.py - repair and apply the sample-type rule in curate_healthy_adults.py

Run it like this (one command, nothing after the file name):

    python scripts/fix_sample_type.py scripts/curate_healthy_adults.py

It makes a backup (.bak), rewrites the sample-type rules with correct
indentation, and checks the file still compiles. Safe to run twice.
"""
import shutil, sys, py_compile

SALIVA = ("saliva|spit|rinse|mouthwash|envo:02000036|oral fluid|"
          "oral cavity|^oral$|^oral |^oral,")
SWAB = r"tongue|swab|buccal|cheek|mouth(?!wash)"
NOTORAL = ("oropharyn|nasopharyn|plaque|stool|skin|water|"
           "neoplasm|carcinoma|cancer|tumou?r")

BLOCK = [
    '    m["sample_type"] = "unknown"',
    f'    m.loc[mat.str.contains("{SALIVA}"), "sample_type"] = "saliva/rinse"',
    f'    m.loc[mat.str.contains(r"{SWAB}"), "sample_type"] = "swab"',
    f'    m.loc[mat.str.contains("{NOTORAL}"), "sample_type"] = "unknown"',
]


def main():
    if len(sys.argv) > 2:
        sys.exit("Too many arguments.\n"
                 "Run exactly:  python scripts/fix_sample_type.py scripts/curate_healthy_adults.py\n"
                 "Then, as a SEPARATE command, run finalize_dataset.sh")
    path = sys.argv[1] if len(sys.argv) == 2 else "scripts/curate_healthy_adults.py"

    try:
        src = open(path).read()
    except FileNotFoundError:
        sys.exit(f"Cannot find {path}. Are you in the ~/oral-age-clock folder?")

    lines = src.split("\n")
    # the block starts at the sample_type initialisation and ends at the last
    # consecutive line that assigns sample_type
    start = None
    for i, l in enumerate(lines):
        if 'm["sample_type"] = "unknown"' in l:
            start = i
            break
    if start is None:
        sys.exit("Could not find the sample_type block in that file. "
                 "Is this the right script?")
    end = start
    for j in range(start + 1, len(lines)):
        s = lines[j].strip()
        if not s:
            break
        if 'sample_type"]' in s or s.startswith(('"', "'")):
            end = j
        else:
            break

    if lines[start:end + 1] == BLOCK:
        print("Already correct - nothing to change.")
        return

    shutil.copy(path, path + ".bak")
    new = lines[:start] + BLOCK + lines[end + 1:]
    open(path, "w").write("\n".join(new))
    try:
        py_compile.compile(path, doraise=True)
    except py_compile.PyCompileError as e:
        sys.exit(f"Still broken:\n{e}\n\nRestore the original with:\n  mv {path}.bak {path}")
    print(f"Fixed. Backup saved as {path}.bak")
    print("The rules are now:")
    for l in BLOCK:
        print("   " + l.strip())
    print("\nNext, run this as a separate command:")
    print("   bash scripts/finalize_dataset.sh metadata/harvest_2026-09-28_1351")


if __name__ == "__main__":
    main()
