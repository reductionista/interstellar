#!/usr/bin/env python3
"""
Create a cleaned observatory code file from ObsCodesNew.txt.

The MPC file uses a fixed-column-ish packed format where numeric fields
run together without delimiters. The C++ read_obscode_file() parser does:
    stream >> code >> obslon >> plxcos >> plxsin
This works for most entries because C++ strtod stops at the second decimal
point and at sign characters. It breaks when:
  1. No numeric coordinates exist (space telescopes).
  2. The observatory name starts with E/e — misread as a scientific-notation
     exponent, corrupting the stream and dropping all subsequent entries.

Fix: replicate strtod semantics in Python to extract the three coordinate
values, then rewrite every line with explicit spaces between all fields.
"""

import os
import re

TESTS   = os.path.join(os.path.dirname(__file__), "../heliolinc2_src/tests")
INFILE  = os.path.join(TESTS, "ObsCodesNew.txt")
OUTFILE = os.path.join(TESTS, "ObsCodes_clean.txt")

# Matches one floating-point number the way C++ strtod would:
# optional sign, digits, optional single decimal point + more digits,
# optional exponent only if followed by a sign/digit (not a letter like E in "Evpatoria").
_FLOAT_RE = re.compile(
    r'[+-]?'            # optional sign
    r'(?:\d+\.?\d*'    # digits with optional decimal
    r'|\.\d+)'          # or leading-decimal form
    r'(?:[eE][+-]?\d+)?' # optional exponent (only if digits follow)
)

def extract_float(s):
    """Return (value, remaining_string) or raise ValueError."""
    s = s.lstrip()
    m = _FLOAT_RE.match(s)
    if not m:
        raise ValueError(f"no float in {s!r:.30}")
    return float(m.group()), s[m.end():]


kept = skipped = c57_remapped = 0

with open(INFILE) as fin, open(OUTFILE, "w") as fout:
    fout.write(fin.readline())  # preserve header
    for line in fin:
        line = line.rstrip("\n")
        if len(line) < 4:
            skipped += 1
            continue

        code = line[0:3]
        rest = line[3:]

        # C57 is listed as TESS but early 3I/ATLAS used it for ATLAS-MLO.
        if code.strip() == "C57":
            fout.write("C57 204.42395 0.943290 +0.332467 ATLAS-MLO, Mauna Loa (old code)\n")
            c57_remapped += 1
            kept += 1
            continue

        try:
            lon, rest = extract_float(rest)
            cos, rest = extract_float(rest)
            sin, rest = extract_float(rest)
        except ValueError:
            skipped += 1  # space telescope or malformed entry
            continue

        name = rest.lstrip()
        # Rewrite with explicit spaces — the C++ >> parser can't mistake
        # a leading 'E' in the name for a scientific-notation exponent.
        fout.write(f"{code} {lon:.6f} {cos:.6f} {sin:+.6f} {name}\n")
        kept += 1

print(f"Kept    {kept} entries")
print(f"Skipped {skipped} entries (space telescopes / malformed)")
print(f"C57 remapped to ATLAS-MLO: {c57_remapped}")
print(f"Output: {OUTFILE}")
