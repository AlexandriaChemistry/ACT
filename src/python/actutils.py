import os

def actdata() -> str:
    actdata = "ACTDATA"
    if not actdata in os.environ:
        sys.exit("Please source ACTRC before running this script. No environment variable %s." % actdata)
    return os.environ[actdata]

def act_library_filename(filename:str) -> str:
    if not os.path.exists(filename):
        filename = ("%s/%s" % ( actdata(), filename ))
    return filename

