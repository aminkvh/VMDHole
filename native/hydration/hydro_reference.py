#!/usr/bin/env python3
"""Reference implementation of the plugin's hydrophobicity colouring.

This is a faithful port of ::VMDHole::colorize_hydrophobic (the built-in Tcl
path) used to validate that sos_triangle_fast's compiled --hydro output is
identical. It also writes the per-frame atom sidecar exactly as the plugin's
write_hydro_sidecar does, so the same file feeds both implementations.

Usage:
  hydro_reference.py sidecar  SPH PDB SCHEME > atoms.dat
  hydro_reference.py colorize SPH ATOMS BASE_PLOT SCHEME > colors_per_triangle.txt

The KD/WW scales and colour thresholds below MUST match vmdhole.tcl.
"""
import sys, math

KD = {'ILE':4.5,'VAL':4.2,'LEU':3.8,'PHE':2.8,'CYS':2.5,'MET':1.9,'ALA':1.8,
      'GLY':-0.4,'THR':-0.7,'SER':-0.8,'TRP':-0.9,'TYR':-1.3,'PRO':-1.6,
      'HIS':-3.2,'HSE':-3.2,'HSD':-3.2,'HSP':-3.2,'HID':-3.2,'HIE':-3.2,'HIP':-3.2,
      'GLN':-3.5,'ASN':-3.5,'GLU':-3.5,'GLUP':-3.5,'ASP':-3.5,'ASPP':-3.5,
      'LYS':-3.9,'ARG':-4.5}
WW = {'TRP':1.85,'PHE':1.13,'TYR':0.94,'LEU':0.56,'ILE':0.31,'MET':0.23,
      'CYS':-0.24,'VAL':0.07,'ALA':-0.17,'GLY':-0.01,'SER':-0.13,'THR':-0.14,
      'PRO':-0.45,'ASN':-0.42,'GLN':-0.58,'HIS':-0.17,'HSE':-0.17,'HSD':-0.17,
      'HSP':-0.17,'HID':-0.17,'HIE':-0.17,'HIP':-0.17,'ARG':-0.81,'LYS':-0.99,
      'ASP':-1.23,'ASPP':-1.23,'GLU':-2.02,'GLUP':-2.02}

def kd(r): return KD.get(r.upper(), 0.0)
def ww(r): return WW.get(r.upper(), 0.0)

def color(h, scheme):
    if scheme == 'kd':
        return ('blue' if h<-3.0 else 'iceblue' if h<-1.5 else 'cyan' if h<-0.5
                else 'white' if h<0.5 else 'yellow' if h<1.5 else 'orange' if h<3.0 else 'red')
    return ('blue' if h<-1.5 else 'iceblue' if h<-0.7 else 'cyan' if h<-0.2
            else 'white' if h<0.2 else 'yellow' if h<0.7 else 'orange' if h<1.2 else 'red')

def read_spheres(path):
    sp = []
    for line in open(path):
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        try:
            sp.append([float(line[30:38]), float(line[38:46]),
                       float(line[46:54]), float(line[60:66])])
        except ValueError:
            pass
    return sp

def bbox(sp):
    ext = max(s[3] for s in sp) + 8.0
    return (min(s[0] for s in sp)-ext, max(s[0] for s in sp)+ext,
            min(s[1] for s in sp)-ext, max(s[1] for s in sp)+ext,
            min(s[2] for s in sp)-ext, max(s[2] for s in sp)+ext)

def r3(v):  # match "%8.3f" the surface is written with
    return round(v, 3)

def cmd_sidecar(sph, pdb, scheme):
    sp = read_spheres(sph)
    xlo,xhi,ylo,yhi,zlo,zhi = bbox(sp)
    for line in open(pdb):
        if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
            continue
        try:
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        if xlo<x<xhi and ylo<y<yhi and zlo<z<zhi:
            res = line[17:21].strip()
            print("%g %g %g %g %g" % (x,y,z,kd(res),ww(res)))

def cmd_colorize(sph, atoms_path, base_plot, scheme):
    sp = read_spheres(sph)
    side = []
    for line in open(atoms_path):
        t = line.split()
        if len(t) >= 5:
            side.append(tuple(float(v) for v in t[:5]))
    col = 3 if scheme == 'kd' else 4
    for s in sp:
        p2 = (s[3]+3.0)**2; tot = 0.0; cnt = 0
        for a in side:
            d2 = (a[0]-s[0])**2 + (a[1]-s[1])**2 + (a[2]-s[2])**2
            if d2 <= p2:
                tot += a[col]; cnt += 1
        s.append(tot/cnt if cnt else 0.0)
    n = len(sp)
    if n > 200:
        stride = math.ceil(n/200.0)
        thin = [sp[i] for i in range(0, n, stride)]
        if thin[-1] is not sp[-1]:
            thin.append(sp[-1])
    else:
        thin = sp
    for line in open(base_plot):
        t = line.split()
        if len(t) >= 2 and t[0] == 'draw' and t[1] == 'trinorm':
            nums = [float(v) for v in t if v not in ('{','}','draw','trinorm')]
            v1, v2, v3 = nums[0:3], nums[3:6], nums[6:9]
            tx = (r3(v1[0])+r3(v2[0])+r3(v3[0]))/3.0
            ty = (r3(v1[1])+r3(v2[1])+r3(v3[1]))/3.0
            tz = (r3(v1[2])+r3(v2[2])+r3(v3[2]))/3.0
            best, best_h = 1e20, 0.0
            for s in thin:
                d = (tx-s[0])**2 + (ty-s[1])**2 + (tz-s[2])**2
                if d < best:
                    best, best_h = d, s[4]
            print(color(best_h, scheme))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'sidecar':
        cmd_sidecar(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == 'colorize':
        cmd_colorize(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        sys.exit("unknown mode: " + mode)
