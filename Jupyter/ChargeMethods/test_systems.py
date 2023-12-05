import numpy as np

def get_system(molname):
    mol = None
    if molname == "carbon-monoxide":
        bonds = [ np.array([(0,1)]), 
                  np.array([(1,0)]) ]
        mol = { "names": np.array([ "c2", "o" ]),
                "coords": np.array([[ 0.0000,    0.0000,    0.060 ],
                                    [ 0.0000,    0.0000,   -0.060 ]]),
                "qtotal": 0,
                "qshell": np.array([ -2, -2 ]),
                "atomnr": np.array([ 6, 8 ]),
                "bonds": bonds }
    elif molname == "water":
        bonds = [ np.array([(0,1),(1,2)]),
                  np.array([(0,1),(2,1)]),
                  np.array([(1,0),(1,2)]) ]
        mol = { "names": np.array([ "hw", "ow", "hw"]),
                "coords": np.array([[ 0.0000,    0.7634,   -0.4681 ],
                                    [ 0.0000,    0.0000,    0.1170 ],
                                    [ 0.0000,   -0.7634,   -0.4681 ]]),
                "qtotal": 0,
                "qshell": np.array([ -1, -2, -1 ]),
                "atomnr": np.array([ 1, 8, 1 ]),
                "bonds": bonds }
    elif molname == "methanol":
        bonds = [ np.array([(0,2),(1,0),(2,4),(3,2),(5,2)]),
                  np.array([(0,1),(0,2),(2,3),(2,4),(2,5)]),
                  np.array([(0,1),(3,2),(0,2),(2,4),(5,2)]),
                  np.array([(2,0),(0,1),(2,3),(2,4),(5,2)]),
                  np.array([(0,2),(1,0),(2,4),(5,2),(3,2)]) ]
        mol = { "names": np.array([ "oh", "hp", "c3", "h1", "h1", "h1" ]),
                "coords": np.array([[  0.7487,    0.1221,    0.0000 ],
                                    [  1.1514,   -0.7502,   -0.0001 ],
                                    [ -0.6672,   -0.0205,    0.0000 ],
                                    [ -1.0269,   -0.5448,   -0.8905 ],
                                    [ -1.0269,   -0.5438,    0.8910 ],
                                    [ -1.0837,    0.9848,   -0.0006 ]]),
                "qtotal": 0,
                "qshell": np.array([ -2, -1, -2, -1, -1, -1 ]),
                "atomnr": np.array([ 8, 1, 6, 1, 1, 1 ]),
                "bonds": bonds}
    elif molname == "acetate":
        bonds = [ np.array([(0,1),(0,4),(0,5),(0,6),(1,2),(1,3)]),
                  np.array([(0,4),(1,0),(0,5),(0,6),(3,1),(1,2)]) ]
        mol = { "names": np.array([ "c3", "cm", "om", "om", "hc", "hc", "hc"]),
                "coords": np.array([[  1.294,  -0.036,   0.384 ],
                                    [ -0.198,   0.010,  -0.068 ],
                                    [ -0.845,   1.021,   0.297 ],
                                    [ -0.596,  -0.981,  -0.727 ],
                                    [  1.702,   0.971,   0.496 ],
                                    [  1.347,  -0.529,   1.360 ],
                                    [  1.899,  -0.615,  -0.316 ]]),
                "qtotal": -1,
                "qshell": np.array([ -2, -2, -2, -2, -1, -1, -1 ]),
                "atomnr": np.array([ 6, 6, 8, 8, 1, 1, 1]),
                "bonds": bonds }
    
    return mol
