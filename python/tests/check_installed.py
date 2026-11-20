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
@brief Check of an installed package act, run by install_test.cmake.

Usage: check_installed.py SITE_DIR PROGRAM, with the directory in which the
package was installed and the path at which the installation puts the
alexandria program. The program itself need not be there: this script puts
a stand-in at that path.
"""
import os
import stat
import sys


def main() -> int:
    """!
    @brief Import the installed package and check that it finds the program
    of its installation.
    @return 0 if all is well, 1 with a message otherwise.
    """
    site_dir, program = (os.path.realpath(name) for name in sys.argv[1:3])
    sys.path.insert(0, site_dir)
    import act  # pylint: disable=import-outside-toplevel

    package_dir = os.path.dirname(os.path.realpath(act.__file__))
    if package_dir != os.path.join(site_dir, "act"):
        print("Imported %s, which is not the installed package in %s"
              % (act.__file__, site_dir))
        return 1
    for name in ("commands", "errors", "legacy", "options", "runner", "_installed"):
        if not os.path.isfile(os.path.join(package_dir, name + ".py")):
            print("The module %s was not installed" % name)
            return 1

    from act import runner  # pylint: disable=import-outside-toplevel

    # Without the program there is nothing to find
    if os.path.exists(program):
        os.remove(program)
    if runner.installed_program() is not None:
        print("Found a program that is not there: %s" % runner.installed_program())
        return 1
    os.makedirs(os.path.dirname(program), exist_ok=True)
    with open(program, "w", encoding="utf-8") as handle:
        handle.write("#!/bin/sh\nexit 0\n")
    os.chmod(program, os.stat(program).st_mode | stat.S_IXUSR)
    found = runner.installed_program()
    if found is None or os.path.realpath(found) != program:
        print("The package looks for its program at %s, the installation puts it at %s"
              % (found, program))
        return 1
    # The program of the installation goes before one in the PATH, the
    # environment variable before both
    os.environ.pop("ACT_ALEXANDRIA", None)
    words = act.config.executable_words()
    if os.path.realpath(words[0]) != program:
        print("The package would run %s instead of %s" % (words[0], program))
        return 1
    os.environ["ACT_ALEXANDRIA"] = sys.executable
    words = act.config.executable_words()
    if os.path.realpath(words[0]) != os.path.realpath(sys.executable):
        print("ACT_ALEXANDRIA does not go before the installed program")
        return 1
    print("The package act in %s finds %s" % (package_dir, program))
    return 0


if __name__ == "__main__":
    sys.exit(main())
