import os
import sys

def actdata() -> str:
    """The directory with the data files of ACT: that in the environment
    variable ACTDATA, if set, and otherwise share/act of the installation or
    share of the source tree that this file is part of."""
    actdata = "ACTDATA"
    if actdata in os.environ and os.environ[actdata]:
        return os.environ[actdata]
    folder = os.path.dirname(os.path.realpath(__file__))
    while True:
        for sub in [ os.path.join("share", "act"), "share" ]:
            datadir = os.path.join(folder, sub)
            if os.path.exists(os.path.join(datadir, "atom_bond.xml")):
                return datadir
        parent = os.path.dirname(folder)
        if parent == folder:
            sys.exit("Cannot find the data files of ACT. Please set the environment variable %s to their directory." % actdata)
        folder = parent

def act_library_filename(filename:str) -> str:
    if not os.path.exists(filename):
        filename = ("%s/%s" % ( actdata(), filename ))
    return filename

