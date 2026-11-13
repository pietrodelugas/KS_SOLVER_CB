#!/usr/bin/env python3
"""max |difference| (eV) between a qeToy eigenvalue file and the eigenvalues of a pw.x data-file-schema.xml.
usage: compare_xml.py qeToy_eigenvalues.txt data-file-schema.xml"""
import re, sys
def xml(path):
    x = open(path).read()
    return [[float(v)*27.211386245988 for v in re.search(r'<eigenvalues[^>]*>(.*?)</eigenvalues>', b, flags=re.S).group(1).split()]
            for b in re.findall(r'<ks_energies>.*?</ks_energies>', x, flags=re.S)]
def txt(path):
    return [[float(v) for v in l.split()[4:]] for l in open(path) if not l.startswith('#') and l.strip()]
a, b = txt(sys.argv[1]), xml(sys.argv[2])
if len(a) == 2*len(b):   # LSDA: the xml keeps both spins in one list per k-point (up bands, then down bands); the file lists all up, then all down
    n = len(b); a = [a[i] + a[i+n] for i in range(n)]
assert len(a) == len(b), (len(a), len(b))
d = [max(abs(x-y) for x, y in zip(r, s)) for r, s in zip(a, b)]
print("%d k-points, max |diff| = %.2e eV (at k-point %d)" % (len(a), max(d), d.index(max(d))+1))
