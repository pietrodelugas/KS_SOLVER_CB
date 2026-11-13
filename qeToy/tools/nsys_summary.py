#!/usr/bin/env python3
"""Summary of an Nsight Systems report of qe_davidson.x / cb_davidson.x: kernel launches, kernel time, GPU busy fraction
(union over all streams) for the whole trace and for its last half (steady state), and the most frequent kernels.
usage: nsys_summary.py report.nsys-rep      (run first: nsys profile -t cuda,nvtx -o report ./qe_davidson.x < input)"""
import sys, csv, subprocess
rep = sys.argv[1]
def stats(name):
    out = subprocess.run(['nsys', 'stats', '--force-export=true', '-r', name, '--format', 'csv', rep], capture_output=True, text=True).stdout
    return [l for l in csv.reader(out.splitlines())]
k = [l for l in stats('cuda_gpu_kern_sum') if l and l[0].replace('.', '').isdigit()]
tot = sum(int(x[1]) for x in k); n = sum(int(x[2]) for x in k)
print('kernel launches %d, kernel time %.3f s, mean kernel %.1f us' % (n, tot/1e9, tot/max(n, 1)/1e3))
k.sort(key=lambda x: -int(x[2]))
for x in k[:6]: print('  %7d launches %6.1f%% of kernel time  %s' % (int(x[2]), 100*int(x[1])/tot, x[-1][:70]))
rows = stats('cuda_gpu_trace')
h = next(i for i, l in enumerate(rows) if l and l[0].startswith('Start'))
ci = {nme: i for i, nme in enumerate(rows[h])}
iv = sorted((float(l[0]), float(l[0]) + float(l[1])) for l in rows[h+1:] if l and l[0].replace('.', '').isdigit())
t0, t1 = iv[0][0], max(e for _, e in iv)
def busy(a):
    cur = a; tot = 0
    for s, e in iv:
        if e <= a: continue
        s = max(s, cur, a)
        if e > s: tot += e - s; cur = e
    return tot
for f in (0.0, 0.5):
    a = t0 + f*(t1-t0)
    print('GPU busy from %3.0f%% of the trace: %.1f%% of %.3f s' % (100*f, 100*busy(a)/(t1-a), (t1-a)/1e9))
