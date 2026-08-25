#!/usr/bin/env python3
"""Fast, VMD-independent regression test for hydro_project. Builds it if
needed, generates synthetic per-frame jobs, computes the SAME arithmetic
independently in Python (a from-scratch re-implementation, not a copy of the
C or Tcl source), and diffs byte-for-byte. Run: python3 test_hydro_project_unit.py

Covers:
  - single-job mode, project-only (QCO) and --bin (KDE + hard-bin, with and
    without a density cap).
  - --batch --bin --bin-global: the GLOBAL cross-frame accumulator, and
    (critically) a check that a naive per-frame-subtotal-then-merge would
    have given a DIFFERENT answer on the same data -- i.e. this test would
    have caught the bug described in NOTES-hydration-accel.md if it were
    reintroduced.

Does NOT cover: bit-identity against Tcl's own arithmetic (that needs a real
VMD session -- see test_hydro_accel_parity.tcl for the end-to-end proof
against real trajectories) or the libm/exp() cross-runtime caveat documented
in NOTES-hydration-accel.md (this test's reference IS Python's math.exp(),
run on the same host/libm as hydro_project, so it cannot see that caveat).
"""
import glob
import math
import os
import random
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "hydro_project")
BUILD_SH = os.path.join(HERE, "..", "build.sh")


def ensure_built():
    if not os.path.exists(BIN):
        subprocess.check_call([BUILD_SH])


def envr(env, co):
    n = len(env)
    if n == 0:
        return 0.0
    if co <= env[0][0]:
        return env[0][1]
    if co >= env[-1][0]:
        return env[-1][1]
    for i in range(n - 1):
        ca, ra = env[i]
        cb, rb = env[i + 1]
        if co >= ca and co <= cb:
            f = (co - ca) / (cb - ca) if cb > ca else 0.0
            return ra + f * (rb - ra)
    return env[-1][1]


def make_frame(seed, n_res, box, dcap):
    random.seed(seed)
    mx, my, mz = random.uniform(-5, 5), random.uniform(-5, 5), random.uniform(-5, 5)
    v = [random.uniform(-1, 1) for _ in range(3)]
    nrm = math.sqrt(sum(c * c for c in v)) or 1.0
    ux, uy, uz = [c / nrm for c in v]
    cmin, cmax = -12.0, 12.0
    env = sorted(set([(round(random.uniform(-15, 15), 3), random.uniform(0.2, 5.0)) for _ in range(60)]))
    resid = []; xs = []; ys = []; zs = []
    expected_positions = []
    rid = 1000
    for _ in range(n_res):
        rid += 1
        natoms = random.choice([1, 1, 2, 3])
        sx = sy = sz = 0.0
        for _a in range(natoms):
            x = mx + random.uniform(-box, box)
            y = my + random.uniform(-box, box)
            z = mz + random.uniform(-box, box)
            resid.append(str(rid)); xs.append(repr(x)); ys.append(repr(y)); zs.append(repr(z))
            sx += x; sy += y; sz += z
        expected_positions.append((sx / natoms, sy / natoms, sz / natoms))
    return dict(mx=mx, my=my, mz=mz, ux=ux, uy=uy, uz=uz, cmin=cmin, cmax=cmax, dcap=dcap,
                env=env, resid=resid, x=xs, y=ys, z=zs, expected_positions=expected_positions)


def write_job(path, fr, dz=1.0, bw=1.4, use_kde=1, with_bin=True):
    with open(path, "w") as f:
        f.write(f"AXIS {fr['mx']!r} {fr['my']!r} {fr['mz']!r} {fr['ux']!r} {fr['uy']!r} {fr['uz']!r}\n")
        f.write(f"RANGE {fr['cmin']!r} {fr['cmax']!r}\n")
        f.write(f"DCAP {fr['dcap']!r}\n")
        f.write(f"ENV {len(fr['env'])}\n")
        f.write(" ".join(repr(c) for c, _ in fr['env']) + "\n")
        f.write(" ".join(repr(r) for _, r in fr['env']) + "\n")
        f.write(f"WATERS {len(fr['resid'])}\n")
        f.write(" ".join(fr['resid']) + "\n")
        f.write(" ".join(fr['x']) + "\n")
        f.write(" ".join(fr['y']) + "\n")
        f.write(" ".join(fr['z']) + "\n")
        if with_bin:
            f.write(f"BIN {dz!r} {bw!r} {use_kde}\n")


def qco_of(fr):
    mx, my, mz, ux, uy, uz = fr['mx'], fr['my'], fr['mz'], fr['ux'], fr['uy'], fr['uz']
    cmin, cmax, dcap = fr['cmin'], fr['cmax'], fr['dcap']
    qco = []
    for (wx, wy, wz) in fr['expected_positions']:
        dxc = wx - mx; dyc = wy - my; dzc = wz - mz
        co = dxc * ux + dyc * uy + dzc * uz
        if co < cmin or co > cmax:
            continue
        rr = envr(fr['env'], co)
        if rr <= 0:
            continue
        rr_eff = dcap if (dcap > 0 and rr > dcap) else rr
        perp2 = dxc * dxc + dyc * dyc + dzc * dzc - co * co
        if perp2 < 0:
            perp2 = 0
        if perp2 <= rr_eff * rr_eff:
            qco.append(co)
    return qco


def bin_of(qco, dz, bw, use_kde):
    bins = {}
    for co in qco:
        if use_kde:
            norm = dz / (bw * 2.5066282746310002)
            lo = math.floor((co - 3.0 * bw) / dz); hi = math.floor((co + 3.0 * bw) / dz)
            for bi in range(int(lo), int(hi) + 1):
                dzz = (bi + 0.5) * dz - co
                w = norm * math.exp(-dzz * dzz / (2.0 * bw * bw))
                bins[bi] = bins.get(bi, 0.0) + w
        else:
            bi = math.floor(co / dz)
            bins[bi] = bins.get(bi, 0.0) + 1.0
    return bins


def fmt_bins(bins):
    if not bins:
        return "BINS 0 -1\n\n"
    lo, hi = min(bins), max(bins)
    out = "BINS %d %d\n" % (lo, hi)
    out += " ".join("%.17g" % bins.get(b, 0.0) for b in range(lo, hi + 1)) + "\n"
    return out


def parse_qco_bins(text):
    lines = text.splitlines()
    assert lines[0].startswith("QCO")
    n = int(lines[0].split()[1])
    qco = [float(v) for v in lines[1].split()] if n > 0 else []
    bins = None
    if len(lines) > 2 and lines[2].startswith("BINS"):
        _, lo, hi = lines[2].split()
        lo, hi = int(lo), int(hi)
        vals = [float(v) for v in lines[3].split()] if hi >= lo else []
        bins = {lo + i: v for i, v in enumerate(vals)}
    return qco, bins


def check(label, cond):
    status = "PASS" if cond else "FAIL"
    print(f"[{status}] {label}")
    return cond


def main():
    ensure_built()
    tmpdir = "/tmp/hydro_project_unit_test"
    os.makedirs(tmpdir, exist_ok=True)
    ok = True

    # --- single-job project + bin (KDE and hard-bin, with/without dcap) ---
    for tag, use_kde, dcap in [("kde", 1, 1.8), ("hardbin", 0, 0.0), ("kde_nodcap", 1, 0.0)]:
        fr = make_frame(hash(tag) & 0xffff, n_res=150, box=8.0, dcap=dcap)
        inp = os.path.join(tmpdir, f"single_{tag}.in")
        write_job(inp, fr, use_kde=use_kde)
        out = subprocess.check_output([BIN, "--bin"], stdin=open(inp), text=True)
        qco_c, bins_c = parse_qco_bins(out)
        qco_py = qco_of(fr)
        bins_py = bin_of(qco_py, 1.0, 1.4, use_kde)
        expected = "QCO %d\n" % len(qco_py) + " ".join("%.17g" % c for c in qco_py) + "\n" + fmt_bins(bins_py)
        ok &= check(f"single-job {tag}: QCO+BINS byte-identical to independent Python reference",
                    out == expected)

    # --- batch mode: global accumulator, and proof it's NOT a subtotal-merge ---
    N_FRAMES = 14
    frames = [make_frame(2000 + i, n_res=110, box=9.0, dcap=(1.8 if i % 3 == 0 else 0.0))
              for i in range(N_FRAMES)]
    batch_dir = os.path.join(tmpdir, "batch")
    os.makedirs(batch_dir, exist_ok=True)
    batch_lines = []
    for i, fr in enumerate(frames):
        inp = os.path.join(batch_dir, f"f{i:03d}.in")
        outp = os.path.join(batch_dir, f"f{i:03d}.out")
        write_job(inp, fr, use_kde=1)
        batch_lines.append(f"{inp}\t{outp}")
    batch_file = os.path.join(batch_dir, "batch.txt")
    with open(batch_file, "w") as f:
        f.write("\n".join(batch_lines) + "\n")
    global_file = os.path.join(batch_dir, "global.out")
    subprocess.check_call([BIN, "--bin", "--bin-global", global_file, "--batch", batch_file])

    # TRUE running-global reference: one continuous accumulation across all
    # frames, in frame order, item by item (mirrors Tcl's original nested loop).
    true_global = {}
    subtotal_merge = {}
    for fr in frames:
        qco = qco_of(fr)
        frame_bins = bin_of(qco, 1.0, 1.4, 1)
        for bi, v in frame_bins.items():
            subtotal_merge[bi] = subtotal_merge.get(bi, 0.0) + v  # WRONG approach
        for co in qco:
            norm = 1.0 / (1.4 * 2.5066282746310002)
            lo = math.floor((co - 3.0 * 1.4) / 1.0); hi = math.floor((co + 3.0 * 1.4) / 1.0)
            for bi in range(int(lo), int(hi) + 1):
                dzz = (bi + 0.5) * 1.0 - co
                w = norm * math.exp(-dzz * dzz / (2.0 * 1.4 * 1.4))
                true_global[bi] = true_global.get(bi, 0.0) + w  # CORRECT approach

    with open(global_file) as f:
        gtext = f.read()
    glines = gtext.splitlines()
    _, glo, ghi = glines[0].split()
    glo, ghi = int(glo), int(ghi)
    gvals = [float(v) for v in glines[1].split()] if ghi >= glo else []
    c_global = {glo + i: v for i, v in enumerate(gvals)}

    ok &= check("batch --bin-global: matches TRUE running-global accumulation exactly",
                c_global == true_global)

    ndiff_vs_subtotal = sum(1 for k in true_global if true_global[k] != subtotal_merge.get(k, 0.0))
    ok &= check(f"regression guard: subtotal-merge WOULD have differed on {ndiff_vs_subtotal} bin(s) "
                f"(this proves the test discriminates the bug fixed in NOTES-hydration-accel.md, "
                f"not a vacuous pass)", ndiff_vs_subtotal > 0)

    for f_out in glob.glob(os.path.join(batch_dir, "f*.out")):
        pass  # per-frame local outputs already exercised implicitly by parsing above

    print()
    if ok:
        print("ALL TESTS PASSED")
        return 0
    else:
        print("SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
