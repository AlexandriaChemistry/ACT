import numpy as np
import scipy
import scipy.optimize as optimize
import math, sys

debug        = False
# Physical constant from https://en.wikipedia.org/wiki/Coulomb_constant
ONE_4PI_EPS0 = 14.3996 # eV / e^2
PI           = 3.1415926

# The radii are from Cordero et al. https://doi.org/10.1039/B801115J
def get_radius():
    return { "hc": 0.31, "h1": 0.31, "ho": 0.31, "hw": 0.31, "c3": 0.76, "cm": 0.76, "oh": 0.66,  "om": 0.66 }

def to_atom(word):
    pp = word.split("#")
    if len(pp) == 1:
        ww = word.split("_")
        if len(ww) == 2:
            return ww[0]
    elif len(pp) == 3:
        return (pp[0].split("_")[0], pp[1].split("_")[0])
    sys.exit("Don't know how to handle " + word)
    
model_data = {}
def read_ff(model):
    filename = model + ".csv"
    with open(filename, "r") as inf:
        for line in inf:
            words = line.strip().split("|")
            if len(words) == 3:
                # The force field files contain
                # type|identifier|value
                atype = to_atom(words[1])
                if not words[0] in model_data:
                    model_data[words[0]] = {}
                value = float(words[2])
                if words[0] == "zeta":
                    value = value*0.1
                model_data[words[0]][atype] = value

def get_beta(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["zeta"]

def get_pol(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["alpha"]

def get_chi(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["chi"]

def get_eta(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["jaa"]

def get_delta_eta(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["hardness"]

def get_delta_chi(model):
    if len(model_data.keys()) == 0:
        read_ff(model)
    return model_data["electronegativity"]

def EFcoulomb(xi, xj, namei, namej, model):
    # Physical constant from https://en.wikipedia.org/wiki/Coulomb_constant
    dx     = xi - xj
    dx2    = np.dot(dx,dx.T)
    rij    = math.sqrt(dx2)
    beta   = get_beta(model)
    betaij = beta[namei]*beta[namej]/np.sqrt(beta[namei]**2 + beta[namej]**2)

    ecoul  = ONE_4PI_EPS0*math.erf(betaij*rij)/rij
    fscal  = ONE_4PI_EPS0*(2*betaij*rij*math.exp(-(betaij*rij)**2)/math.sqrt(PI) - math.erf(betaij*rij))/dx2
    return ecoul, (-fscal/rij)*dx
            
def calcJEEM(names, coords, model):
    N = coords.shape[0]
    J = np.zeros((N+1,N+1), dtype=float)
    # Half a matrix of Coulomb interactions
    for i in range(N):
        for j in range(i+1,N):
            E, F = EFcoulomb(coords[i], coords[j], names[i], names[j], model)
            J[i,j] = J[j,i] = 0.5*E
    # Fill the diagonal
    eta  = get_eta(model)
    for i in range(N):
        if names[i] in eta:
            J[i,i] = eta[names[i]]
        else:
            print("No eta for %s" % names[i])
            exit(1)
    # Extra row for enforcing total charge is equal to whatever it should be
    for i in range(N):
        J[N,i] = 1
    # Extra columne for enforcing equal electronegativity
    for i in range(N):
        J[i,N] = -1
    return J

def printDipole(q, coords, atomnumber, have_shells, qshell=None, shellX=None):
    # Determine center of charge
    coq    = np.zeros(3)
    atntot = 0.0
    for i in range(coords.shape[0]):
        coq    += atomnumber[i]*coords[i]
        atntot += atomnumber[i]
    coq = coq/atntot
    DEBYE = 4.80321
    mu    = np.zeros(3)
    for i in range(coords.shape[0]):
        mu   += DEBYE*q[i]*(coords[i]-coq)
        if have_shells:
            mu   += DEBYE*qshell[i]*(shellX[i]-coq)
    mutot = math.sqrt(np.dot(mu,mu.T))
    print("Dipole = %g Debye" % (mutot))

def calcJcs(i, names, coords, shellX, qshell, model):
    N        = coords.shape[0]
    MyshellX = np.reshape(np.array(shellX), (-1, 3))
    if MyshellX.shape[0] < N or names.shape[0] < N or qshell.shape[0] < N:
        print(i)
        print("names {}".format(names))
        print("coords {}".format(coords))
        print("MyshellX {}".format(MyshellX))
        print("qshell {}".format(qshell))
        sys.exit("stop")
    Jcs      = 0
    for j in range(N):
        if i != j:
            EC, F = EFcoulomb(coords[i], MyshellX[j], names[i], names[j], model)
            Jcs += qshell[j]*EC
    return Jcs

def get_bondhardness(ai, aj, model):
    bij = (ai,aj)
    bji = (aj,ai)
    bcc = get_delta_eta(model)
    if bij in bcc:
        return bcc[bij]
    elif bji in bcc:
        return bcc[bji]
    else:
        print("No hardness information for bond %s-%s" % (ai, aj))
        return 0

def getDeltaChi(ai, aj, model):
    bij = (ai,aj)
    bji = (aj,ai)
    bcc = get_delta_chi(model)
    if bij in bcc:
        return bcc[bij]
    elif bji in bcc:
        return -bcc[bji]
    else:
        print("No deltaChi information for bond %s-%s" % (ai, aj))
        return 0

def calcJSQE(names, coords, bonds, verbose, model):
    # Abuse the EEM routine for some calculations
    JEEM  = calcJEEM(names, coords, model)
    Natom = coords.shape[0]
    Nbond = len(bonds)
    if verbose:
        print("There are %d atoms and %d bonds" % (Natom, Nbond))
    JSQE = np.zeros((Nbond,Nbond), dtype=float)
    for n in range(Nbond):
        bi = bonds[n][0]
        bj = bonds[n][1]
        for p in range(Nbond):
            bk = bonds[p][0]
            bl = bonds[p][1]
            JSQE[n][p] = (JEEM[bi][bk] - JEEM[bi][bl] - JEEM[bj][bk] + JEEM[bj][bl])
            if n == p:
                # NOTE: Here bond hardness needs to be added to the atomic hardness that are the 
                # diagonals of JEEM. 
                hardness = get_bondhardness(names[bi], names[bj], model)
                JSQE[n][p] += hardness
                if verbose:
                    print("hardness %g" % hardness)
    if verbose:
        print("JEEM:\n{}".format(JEEM))
        print("JSQE:\n{}".format(JSQE))
    return JEEM, JSQE

def gradient(params, args):
    # Code not correct yet! Use at your own risk.
    N        = args.names.shape[0]
    MyshellX = np.reshape(np.array(params), (-1, 3))
    force    = np.zeros((N,3))
    wall     = 16*args.delta**2
    Etot     = 0
    for i in range(N):
        for j in range(N):
            if j > i:
                Ec_ss, Fc_ss = EFcoulomb(MyshellX[i], MyshellX[j], args.names[i], args.names[j], args.model)
                qq_ss = args.qshell[i]*args.qshell[j]
                force[i] -= Fc_ss*qq_ss
                force[j] += Fc_ss*qq_ss
                Etot     += Ec_ss*qq_ss

            if j != i:
                Ec_qs, Fc_qs = EFcoulomb(args.coords[i], MyshellX[j], args.names[i], args.names[j], args.model)
                qq_qs = args.q[i]*args.qshell[j]
                # We only want the force on the shells!
                force[j] += Fc_qs*qq_qs
                Etot     += Ec_qs*qq_qs

        dx    = MyshellX[i]-args.coords[i]
        kpol  = ONE_4PI_EPS0*(args.qshell[i]**2)/args.pol[args.names[i]]
        force[i] += kpol*dx
        Etot += 0.5*kpol*np.dot(dx, dx.T)
    fflat = np.ndarray.flatten(force)
    if debug:
#        print("shelx {}".format(params))
#        print("force {}".format(np.ndarray.flatten(force)))
        print("RMS force %g" % (np.sqrt(np.dot(fflat, fflat.T))))
    return Etot, fflat

class ShellMinimize:
    def __init__(self, names, coords, qshell, model):
        self.coords = coords
        self.names  = names
        self.pol    = get_pol(model)
        self.qshell = qshell
        self.model  = model
        self.delta  = 0.02 # Angstrom
        N           = self.names.shape[0]
        self.shellX = np.zeros(N*3)
        lb          = np.zeros(N*3)
        ub          = np.zeros(N*3)
        for i in range(N):
            for m in range(3):
                ind = 3*i+m
                xxx = self.coords[i][m]
                self.shellX[ind] = xxx
                lb[ind] = xxx-self.delta
                ub[ind] = xxx+self.delta
        self.bounds = optimize.Bounds(lb, ub, keep_feasible=True)
  
    def minimize(self, q):
        # To be implemented
        self.q = np.copy(q)
        N      = self.names.shape[0]
        result = optimize.minimize(gradient, self.shellX.tolist(), args=self, method="BFGS", jac=True)
        if result.success:
            xdiff     = self.shellX - result.x
            self.rmsd = np.sqrt(np.sum(np.dot(xdiff, xdiff.T))/len(self.names))
            self.shellX = np.copy(result.x)
            return result.fun
        else:
            raise ValueError(result.message)
                
