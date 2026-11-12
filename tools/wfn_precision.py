"""Electron and pair counts implied by the MO coefficients stored in an AIMPAC .wfn file
(single-determinant, closed- or open-shell occupations), computed in float64.

With the 1-RDM D = sum_i n_i c_i c_i^T over the primitive basis and the exact primitive overlap S:
    N'      = tr(DS)
    Coulomb = (tr DS)^2                 (int I^C dR)
    XC      = -1/2 tr(DSDS)  (RHF, n_i = 2) or -tr(D_a S D_a S) - tr(D_b S D_b S) (spin-resolved)
If the MOs were exactly orthonormal these would be N, N^2 and -N. Any deviation is carried by
the input file, not by the intracule code."""
import sys, re, numpy as np
from math import pi, sqrt, exp

NLM = {1:(0,0,0),2:(1,0,0),3:(0,1,0),4:(0,0,1),5:(2,0,0),6:(0,2,0),7:(0,0,2),8:(1,1,0),9:(1,0,1),
       10:(0,1,1),11:(3,0,0),12:(0,3,0),13:(0,0,3),14:(2,1,0),15:(2,0,1),16:(0,2,1),17:(1,2,0),
       18:(1,0,2),19:(0,1,2),20:(1,1,1)}

def fnum(t): return float(t.replace('D','E').replace('d','e'))

def read_wfn(fn):
    L = open(fn).read().split('\n')
    m = re.search(r'(\d+)\s+MOL ORBITALS\s+(\d+)\s+PRIMITIVES\s+(\d+)\s+NUCLEI', L[1])
    nmo, nprim, ncen = map(int, m.groups())
    xyz = []
    for l in L[2:2+ncen]:
        nums = re.findall(r'-?\d+\.\d+', l.split(')')[1])
        xyz.append([float(x) for x in nums[:3]])
    i = 2+ncen
    cen, typ, ex = [], [], []
    while len(cen) < nprim:
        l = L[i]; i += 1; f = l[20:]
        cen += [int(f[k:k+3]) for k in range(0, len(f.rstrip()), 3)]
    while len(typ) < nprim:
        l = L[i]; i += 1; f = l[20:]
        typ += [int(f[k:k+3]) for k in range(0, len(f.rstrip()), 3)]
    while len(ex) < nprim:
        l = L[i]; i += 1; f = l[10:]
        ex += [fnum(f[k:k+14]) for k in range(0, len(f.rstrip()), 14)]
    occ, C = [], []
    while len(occ) < nmo:
        l = L[i]; i += 1
        if not l.startswith('MO'): continue
        occ.append(float(re.search(r'OCC NO =\s*(-?[\d.]+)', l).group(1)))
        c = []
        while len(c) < nprim:
            f = L[i]; i += 1
            c += [fnum(f[k:k+16]) for k in range(0, len(f.rstrip()), 16)]
        C.append(c)
    return np.array(xyz), np.array(cen)-1, np.array(typ), np.array(ex), np.array(occ), np.array(C)

def s1d(i, j, a, b, A, B):
    """1D overlap of (x-A)^i e^{-a(x-A)^2} and (x-B)^j e^{-b(x-B)^2} by Hermite E_0^{ij}."""
    p = a+b; mu = a*b/p; P = (a*A+b*B)/p
    PA, PB = P-A, P-B
    E = {(0,0,0): exp(-mu*(A-B)**2)}
    def e(t, ii, jj):
        if t < 0 or t > ii+jj or ii < 0 or jj < 0: return 0.0
        if (t,ii,jj) in E: return E[(t,ii,jj)]
        if ii > 0:
            v = e(t-1,ii-1,jj)/(2*p) + PA*e(t,ii-1,jj) + (t+1)*e(t+1,ii-1,jj)
        else:
            v = e(t-1,ii,jj-1)/(2*p) + PB*e(t,ii,jj-1) + (t+1)*e(t+1,ii,jj-1)
        E[(t,ii,jj)] = v; return v
    return e(0,i,j)*sqrt(pi/p)

def overlap(xyz, cen, typ, ex):
    n = len(ex); S = np.empty((n,n))
    for p in range(n):
        for q in range(p+1):
            v = 1.0
            for k in range(3):
                v *= s1d(NLM[typ[p]][k], NLM[typ[q]][k], ex[p], ex[q], xyz[cen[p],k], xyz[cen[q],k])
            S[p,q] = S[q,p] = v
    return S

for fn in sys.argv[1:]:
    xyz, cen, typ, ex, occ, C = read_wfn(fn)
    S = overlap(xyz, cen, typ, ex)
    Smo = C @ S @ C.T
    N = round(occ.sum())
    D = (C.T*occ) @ C
    Np = np.trace(D@S)
    coul = Np**2
    if np.allclose(occ, 2.0):
        xc = -0.5*np.trace(D@S@D@S)
    else:
        raise SystemExit('open shell: not implemented')
    tot = coul + xc
    off = np.abs(Smo - np.diag(np.diag(Smo))).max()
    print(f"{fn}\n  N = {N}, MOs = {len(occ)}, primitives = {len(ex)}")
    print(f"  max |<i|i>-1| = {np.abs(np.diag(Smo)-1).max():.3e}   max |<i|j>|, i!=j = {off:.3e}")
    print(f"  N'          = {Np:.10f}   (rel {Np/N-1:+.3e})")
    print(f"  N(N-1)'     = {tot:.8f}   (exact {N*(N-1)})")
    print(f"  N^2'        = {coul:.8f}   (exact {N*N})")
    print(f"  XC'         = {xc:.8f}   (exact {-N})")
