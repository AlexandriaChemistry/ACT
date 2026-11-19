#!/usr/bin/env python3
#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
"""!
@file
@brief A stand-in for the alexandria program in the tests of the Python
interface, such that these tests need no build of ACT.

It knows two commands. `help -export json` writes the description of the
options, as alexandria does. `work` writes its arguments to the JSON file
and behaves as its options say: it can fail, write to standard error or
take a long time.
"""
import json
import os
import sys
import time

## The description of the options that this program gives of itself
DESCRIPTION = {
    "program": "alexandria",
    "format": 1,
    "commands": {
        "work": {
            "description": "Pretend to do a calculation.",
            "help": "alexandria work pretends to do a calculation.\n\nSecond paragraph.",
            "options": [
                {"name": "ff", "type": "file", "description": "Force field", "hidden": False,
                 "required": True, "multiple": False, "default": "aff.xml",
                 "file": {"input": True, "output": False, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".xml"]}},
                {"name": "many", "type": "file", "description": "Several files", "hidden": False,
                 "required": False, "multiple": True, "default": ["in.xml"],
                 "file": {"input": True, "output": False, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".xml"]}},
                {"name": "lib", "type": "file", "description": "Library file", "hidden": False,
                 "required": False, "multiple": False, "default": "lib.csv",
                 "file": {"input": True, "output": False, "library": True,
                          "allow_missing": False, "directory": False, "extensions": [".csv"]}},
                {"name": "c", "type": "file", "description": "Structure", "hidden": False,
                 "required": True, "multiple": False, "default": "confout.pdb",
                 "file": {"input": False, "output": True, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".pdb"]}},
                {"name": "g", "type": "file", "description": "Log", "hidden": False,
                 "required": True, "multiple": False, "default": "work.log",
                 "file": {"input": False, "output": True, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".log"]}},
                {"name": "extra", "type": "file", "description": "Optional output",
                 "hidden": False, "required": False, "multiple": False, "default": "extra.xvg",
                 "file": {"input": False, "output": True, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".xvg"]}},
                {"name": "json", "type": "file", "description": "Results", "hidden": False,
                 "required": False, "multiple": False, "default": "work.json",
                 "file": {"input": False, "output": True, "library": False,
                          "allow_missing": False, "directory": False, "extensions": [".json"]}},
                {"name": "f", "type": "string", "description": "Molecule file", "hidden": False,
                 "required": False, "multiple": False, "default": None},
                {"name": "name", "type": "string", "description": "A name", "hidden": False,
                 "required": False, "multiple": False, "default": None},
                {"name": "maxiter", "type": "int", "description": "Number of steps",
                 "hidden": False, "required": False, "multiple": False, "default": 100},
                {"name": "toler", "type": "real", "description": "Tolerance", "hidden": False,
                 "required": False, "multiple": False, "default": 1e-06},
                {"name": "box", "type": "real", "description": "Box", "hidden": False,
                 "required": False, "multiple": True, "max_values": 3, "default": [0, 0, 0]},
                {"name": "minalg", "type": "enum", "description": "Algorithm", "hidden": False,
                 "required": False, "multiple": False, "default": "LBFGS",
                 "values": ["LBFGS", "Newton", "Steep"]},
                {"name": "oneH", "type": "bool", "description": "One hydrogen type",
                 "hidden": False, "required": False, "multiple": False, "default": False},
                {"name": "del", "type": "bool", "description": "A Python keyword",
                 "hidden": False, "required": False, "multiple": False, "default": False},
                {"name": "fail", "type": "int", "description": "Exit status to fail with",
                 "hidden": True, "required": False, "multiple": False, "default": 0},
                {"name": "sleep", "type": "real", "description": "Seconds to sleep",
                 "hidden": True, "required": False, "multiple": False, "default": 0},
            ],
        },
    },
}


def main(argv):
    """!
    @brief Run the stand-in.
    @param argv The command line, without the name of the program.
    @return The exit status.
    """
    if argv[:3] == ["help", "-export", "json"]:
        with open("alexandria-options.json", "w", encoding="utf-8") as handle:
            json.dump(DESCRIPTION, handle)
        print("rank: 0/1 nodetype: Master")
        return 0
    if not argv or argv[0] != "work":
        sys.stderr.write("'%s' is not an ACT command.\n" % (argv[0] if argv else ""))
        return 1
    # Collect the options: a flag is followed by its values
    options = {}
    flag = None
    for word in argv[1:]:
        is_number = True
        try:
            float(word)
        except ValueError:
            is_number = False
        if word.startswith("-") and not is_number:
            flag = word[1:]
            options[flag] = []
        elif flag is not None:
            options[flag].append(word)
    seconds = float(options.get("sleep", ["0"])[0])
    if seconds > 0:
        time.sleep(seconds)
    status = int(options.get("fail", ["0"])[0])
    log = options.get("g", ["work.log"])[0]
    with open(log, "w", encoding="utf-8") as handle:
        handle.write("log of work\n")
    if status != 0:
        sys.stderr.write("\n-------------------------------------------------------\n"
                         "Program:     alexandria work\n"
                         "Source file: fake.cpp (line 1)\n\n"
                         "Fatal error:\n"
                         "The calculation went wrong.\n"
                         "For more information and tips for troubleshooting, please check\n"
                         "-------------------------------------------------------\n")
        return status
    with open(options.get("c", ["confout.pdb"])[0], "w", encoding="utf-8") as handle:
        handle.write("ATOM\n")
    if "extra" in options:
        with open(options["extra"][0] if options["extra"] else "extra.xvg", "w",
                  encoding="utf-8") as handle:
            handle.write("1 2\n")
    if "json" in options:
        result = {"work": {"argv": argv, "cwd": os.getcwd(),
                           "env": os.environ.get("ACT_TEST_VARIABLE", ""),
                           "rmsd": {"value": 0.0192, "unit": "nm"}}}
        with open(options["json"][0] if options["json"] else "work.json", "w",
                  encoding="utf-8") as handle:
            json.dump(result, handle)
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
