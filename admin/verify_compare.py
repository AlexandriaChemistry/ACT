#!/usr/bin/env python3
"""Compare the regression set of admin/verify.sh with a reference.

For each case the trained force field (the values of all parameters in
ff_inter.xml) and the second virial (all numbers in b2.xvg) are compared.
The largest relative difference of each case is printed, and the exit
status is nonzero if it exceeds the tolerance in any case, or if a file
is missing.
"""

import argparse
import os
import re
import sys

VALUE = re.compile(r'\svalue="([^"]+)"')


def numbers(path):
    """The numbers in a force field or xvg file, in order."""
    values = []
    with open(path) as f:
        if path.endswith(".xml"):
            for line in f:
                if "<parameter " in line:
                    for v in VALUE.findall(line):
                        try:
                            values.append(float(v))
                        except ValueError:
                            pass
        else:
            for line in f:
                if line.startswith(("#", "@")):
                    continue
                values += [float(w) for w in line.split()]
    return values


def relativeDifference(ref, new):
    """The largest relative difference, or None if the numbers do not match up."""
    if len(ref) != len(new):
        return None
    largest = 0.0
    for a, b in zip(ref, new):
        if a != b:
            largest = max(largest, abs(a - b) / max(abs(a), abs(b)))
    return largest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    parser.add_argument("output")
    parser.add_argument("--tolerance", type=float, default=0.0)
    args = parser.parse_args()
    if not os.path.isdir(args.reference):
        print(f"    No reference in {args.reference}, run with --save-reference first")
        return 1
    ok = True
    for case in sorted(os.listdir(args.reference)):
        for name in ("ff_inter.xml", "b2.xvg"):
            ref = os.path.join(args.reference, case, name)
            if not os.path.exists(ref):
                continue
            new = os.path.join(args.output, case, name)
            if not os.path.exists(new):
                print(f"    {case:16s} {name:14s} missing")
                ok = False
                continue
            diff = relativeDifference(numbers(ref), numbers(new))
            if diff is None:
                print(f"    {case:16s} {name:14s} different number of values")
                ok = False
            else:
                status = "identical" if diff == 0 else f"relative difference {diff:.3g}"
                if diff > args.tolerance:
                    status += ", FAILED"
                    ok = False
                print(f"    {case:16s} {name:14s} {status}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
