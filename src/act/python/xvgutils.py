
debugXvgUtils = False
        
def interpret_legend(line:str):
    legval = None

    #find and save subtitle if one exists
    labkey = "subtitle"
    if line.find(labkey) >= 0:
        legval = line[line.find(labkey)+len(labkey)+1:].strip()
        legval = legval[1:-1]
        return labkey, legval

    labkey = "title"
    if line.find(labkey) >= 0:
        legval = line[line.find(labkey)+len(labkey)+1:].strip()
        legval = legval[1:-1]
        return labkey, legval
    
    for axis in [ "x", "y" ]:
        labkey  = axis+"axis"
        labkey2 = " label" #Space is added to avoid ticklegends for old Grace purposes
        if line.find(labkey) >= 0 and line.find(labkey2) >= 0:
            legval = line[line.find(labkey2)+len(labkey2)+1:].strip()
            legval = legval[1:-1]
            return axis+"label", legval
    legkey = "legend"
    if line.find(legkey) >= 0 and line[0] == 's':
        labval = line[line.find(legkey)+len(legkey)+1:].strip()
        labval = labval[1:-1]
        return legkey, labval
    return None, None

class xvgDataSet:
    '''A simple class to hold an xvg data set'''
    def __init__(self, line:str):
        self.x = []
        self.y = []
        self.dy = None
        if line.find("xydy") >= 0:
            self.dy = []
        self.xmin = 1e9
        self.xmax = -1e9
        self.ymin = 1e9
        self.ymax = -1e9

    def have_dy(self)->bool:
        return self.dy != None

    def set_label(self, label:str):
        self.label = label

    def add_point(self, x:float, y:float, dy=None):
        self.x.append(x)
        self.y.append(y)
        self.xmin = min(self.xmin, x)
        self.xmax = max(self.xmax, x)
        if dy:
            self.dy.append(dy)
            self.ymin = min(self.ymin, y-dy)
            self.ymax = max(self.ymax, y+dy)
        else:
            self.ymin = min(self.ymin, y)
            self.ymax = max(self.ymax, y)

def read_xvg(filename:str, residual:bool=False, filelabel:bool=False):
    labels  = {}
    legends   = []
    dataset  = []
    numwords = None
    newset   = None
    with open(filename, "r") as inf:
        for line in inf:
            nhash = line.find("#")
            line = line[:nhash]
            if len(line) == 0:
                continue
            
            nleg = line.find("@")
            if nleg >= 0:
                myline = line[nleg+1:].strip()
                if myline.find("type") >= 0 and myline.find("loctype") < 0:
                    dataset.append(xvgDataSet(myline))
                elif len(myline) > 0:
                    labkey, legval = interpret_legend(myline)
                    if labkey and legval:
                        if labkey == "legend":
                            if filelabel:
                                legval += " " + filename
                            legends.append(legval)
                        else:
                            labels[labkey] = legval
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
                        dataset.append(xvgDataSet(""))
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
                    if len(dataset) > 0 and dataset[0].have_dy() and numwords == 3:
                        try:
                            xx = float(w[0])
                            yy = float(w[1])
                            dy = float(w[2])
                            if residual:
                                yy -= xx
                            dataset[0].add_point(xx, yy, dy)
                        except:
                            if debugXvgUtils:
                                print("Could not read line '%s'" % line)
                    else:
                        if len(dataset) < numwords-1:
                            for i in range(len(dataset), numwords-1):
                                dataset.append(xvgDataSet(""))
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
        if not ylabel in labels:
            labels[ylabel] = ""
        labels[ylabel] += " (Residual)"
    return legends, labels, dataset
    
