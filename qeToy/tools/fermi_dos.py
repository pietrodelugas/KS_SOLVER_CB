#!/usr/bin/env python3
"""Fermi energy and DOS from an eigenvalue file (qeToy or reference layout: kx ky kz weight e1..en, eV).
Marzari-Vanderbilt cold smearing as in QE (wgauss/w0gauss, n=-1).
usage: fermi_dos.py file nelec degauss_Ry [dos_out [other_file]]    (other_file: prints the max DOS and Ef difference)"""
import sys, math
RY = 13.605693122994
def load(p):
    k = []; 
    for l in open(p):
        if l.startswith('#') or not l.strip(): continue
        v = l.split(); k.append((float(v[3]), [float(x) for x in v[4:]]))
    s = sum(w for w, _ in k)
    return [(2.0*w/s, e) for w, e in k]          # weights summing to 2 (spin unpolarised)
def wgauss(x):
    xp = x - 1/math.sqrt(2); return 0.5*(1+math.erf(xp)) + math.exp(-min(200, xp*xp))/math.sqrt(2*math.pi)
def w0gauss(x):
    xp = x - 1/math.sqrt(2); return math.exp(-min(200, xp*xp))/math.sqrt(math.pi)*(2 - math.sqrt(2)*x)
def fermi(k, nel, dg):
    f = lambda ef: sum(w*wgauss((ef-e)/dg) for w, es in k for e in es) - nel
    lo, hi = min(e for _, es in k for e in es), max(e for _, es in k for e in es)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    return 0.5*(lo+hi)
def dos(k, dg, grid): return [sum(w*w0gauss((x-e)/dg) for w, es in k for e in es)/dg for x in grid]
if __name__ == '__main__':
    p, nel, dgry = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]); dg = dgry*RY
    k = load(p); ef = fermi(k, nel, dg)
    grid = [-5 + 0.05*i for i in range(int((45+5)/0.05)+1)]; d = dos(k, dg, grid)
    print("%s: %d k-points, Fermi energy %.4f eV" % (p, len(k), ef))
    if len(sys.argv) > 4:
        with open(sys.argv[4], 'w') as fo:
            fo.write("# E(eV)  DOS(states/eV, both spins)\n")
            for x, y in zip(grid, d): fo.write("%8.3f %12.6f\n" % (x, y))
    if len(sys.argv) > 5:
        k2 = load(sys.argv[5]); ef2 = fermi(k2, nel, dg); d2 = dos(k2, dg, grid)
        print("other: Fermi %.4f eV; |dEf| = %.2e eV, max|dDOS| = %.2e (max DOS %.3f)" % (ef2, abs(ef-ef2), max(abs(a-b) for a, b in zip(d, d2)), max(d2)))
