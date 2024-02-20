
debugXvgUtils = False
        
def interpret_legend(line):
    legval = None
    legkey = "title"
    if line.find(legkey) >= 0:
        legval = line[line.find(legkey)+len(legkey)+1:].strip()
        legval = legval[1:-1]
        return legkey, legval
    
    for axis in [ "x", "y" ]:
        legkey  = axis+"axis"
        legkey2 = "label" 
        if line.find(legkey) >= 0 and line.find(legkey2) >= 0:
            legval = line[line.find(legkey2)+len(legkey2)+1:].strip()
            legval = legval[1:-1]
            return axis+"label", legval
    legkey = "legend"
    if line.find(legkey) >= 0 and line[0] == 's':
        legval = line[line.find(legkey)+len(legkey)+1:].strip()
        legval = legval[1:-1]
        return "label", legval
    return None, None

class xvgDataSet:
    '''A simple class to hold an xvg data set'''
    def __init__(self):
        self.x = []
        self.y = []
        self.xmin = 1e9
        self.xmax = -1e9
        self.ymin = 1e9
        self.ymax = -1e9
        
    def set_label(self, label:str):
        self.label = label

    def add_point(self, x:float, y:float):
        self.x.append(x)
        self.y.append(y)
        self.xmin = min(self.xmin, x)
        self.ymin = min(self.ymin, y)
        self.xmax = max(self.xmax, x)
        self.ymax = max(self.ymax, y)

def read_xvg(filename:str, residual:bool=False, filelabel:bool=False):
    legends  = {}
    labels   = []
    dataset  = []
    numwords = None
    newset   = None
    with open(filename, "r") as inf:
        for line in inf:
            nhash = line.find("#")
            if nhash == 0:
                line = line[:nhash]
                continue
            
            nleg = line.find("@")
            if nleg >= 0:
                myline = line[nleg+1:].strip()
                if line.find("@type") == 0:
                    dataset.append(xvgDataSet())
                elif len(myline) > 0:
                    legkey, legval = interpret_legend(myline)
                    if legkey and legval:
                        if legkey == "label":
                            if filelabel:
                                legval += " " + filename
                            labels.append(legval)
                        else:
                            legends[legkey] = legval
                continue
            
            w = line.split()
            if len(w) == 1:
                w = line.split(",")
            if None == numwords:
                numwords = len(w)
            if len(w) == numwords:
                if numwords == 2:
                    if len(dataset) == 0:
                        if debugXvgUtils:
                            print("Found data but no dataset yet")
                        dataset.append(xvgDataSet())
                    try:
                        xx = float(w[0])
                        yy = float(w[1])
                        if residual:
                            yy -= xx
                        dataset[len(dataset)-1].add_point(xx, yy)
                    except:
                        if debugXvgUtils:
                            print("Could not read line '%s'" % line)
                elif numwords > 2:
                    if len(dataset) < numwords-1:
                        for i in range(len(dataset), numwords-1):
                            dataset.append(xvgDataSet())
                    for i in range(numwords-1):
                        try:
                            xx = float(w[0])
                            yy = float(w[i+1])
                            if residual:
                                yy -= xx
                            dataset[i].add_point(xx, yy)
                        except:
                            if debugXvgUtils:
                                print("Could not read line '%s'" % line)
    if residual:
        ylabel = "ylabel"
        if not ylabel in legends:
            legends[ylabel] = ""
        legends[ylabel] += " (Residual)"
    return labels, legends, dataset
    
