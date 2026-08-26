#!/usr/bin/env python3
"""
Apply the VMDHole parallel-CONNOLLY / fast-CAPSULE patches to a HOLE2 src tree.

Idempotent and non-destructive: every file it changes in place is first backed
up to <file>.vmdhole_orig (only once - the pristine copy is preserved). Running
it twice is a no-op. New source files are simply copied in.

Usage:  apply_patches.py  <hole2/src dir>

What it does (see connolly-parallel-twopass-foundation memory / the .f headers):
  * copies 8 new source files: hcapen_fast.f, coarea_fast.f, holcal_par.f,
    holeen_par.f, sphqpu_par.f, h2dmap_par.f, tsatr_fast.f, concal_par.f
  * machine_dep.g77 : adds an inert RNG draw counter (COMMON /RNGCNT/ NDRAW)
  * Makefile        : swaps the object files for the fast/parallel versions,
                      defines OMPFLAGS, links with -fopenmp, and adds per-file
                      -fopenmp rules for ONLY the parallel-region files
                      (global -fopenmp would push hole.f's 24MB arrays onto the
                      stack and seg-fault every mode).
concal.f is not edited in place: concal_par.f is compiled instead of it (see
OBJ_SWAPS), and differs from upstream in one gated WRITE. Like the other
parallel-region files it is compiled with -fopenmp so its arrays become
per-thread.
"""
import sys, os, re, shutil, hashlib

# sha256 (first 16 hex chars) of the PRISTINE upstream files these patches were
# written against, at HOLE2 a8eaf6121ba66625446933f4acd7d6aa336dbb47. The
# patches REPLACE whole files rather than applying diffs, so a drifted upstream
# would be overwritten silently and the result would be a binary built from a
# mixture nobody validated. Checked, not assumed. Set VMDHOLE_SKIP_BASE_CHECK=1
# to override deliberately.
BASE_SHA16 = {
    "hcapen.f":        "0fbfb016033a12ba",
    "tsatr.f":         "96ed23294b1bc460",
    "coarea.f":        "53747f0e4dd18860",
    "holcal.f":        "9a1f2696dca2a9a2",
    "holeen.f":        "def08fb82ffac6f2",
    "h2dmap.f":        "d6f0e2f40af51377",
    "machine_dep.g77": "e793bb78f9241bf0",
    # Not edited in place - replaced at link time via OBJ_SWAPS. Pinned for the
    # same reason as the rest: a drifted upstream would be shadowed by a stale
    # _par copy and the mixture would never be noticed.
    "concal.f":        "ca0e9cdf400ecf6e",
    "sphqpu.f":        "19341d79302c9440",
    "Makefile":        "40206adfb8852861",
}


def _sha16(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def check_base_files(src):
    """Verify the upstream files we are about to replace are the ones these
    patches were validated against. Compares the .vmdhole_orig backup when one
    exists, so re-running on an already-patched tree still checks the original."""
    if os.environ.get("VMDHOLE_SKIP_BASE_CHECK"):
        print("  base-file check : SKIPPED (VMDHOLE_SKIP_BASE_CHECK set)")
        return
    bad = []
    for f, want in BASE_SHA16.items():
        path = os.path.join(src, f)
        orig = path + ".vmdhole_orig"
        probe = orig if os.path.isfile(orig) else path
        if not os.path.isfile(probe):
            bad.append("%s: missing" % f)
            continue
        got = _sha16(probe)
        if got != want:
            bad.append("%s: got %s, expected %s" % (f, got, want))
    if bad:
        sys.exit(
            "ERROR: upstream sources differ from the revision these patches were\n"
            "       written and validated against:\n         "
            + "\n         ".join(bad)
            + "\n       Build with the pinned revision (build-vmdhole-optimized.sh does\n"
              "       this automatically), or re-validate and update BASE_SHA16 in\n"
              "       this file. Set VMDHOLE_SKIP_BASE_CHECK=1 to override.")
    print("  base-file check : %d file(s) match the pinned upstream" % len(BASE_SHA16))

NEW_FILES = ["hcapen_fast.f", "coarea_fast.f", "holcal_par.f", "holeen_par.f",
             "sphqpu_par.f", "h2dmap_par.f", "tsatr_fast.f", "concal_par.f"]
# object-file swaps in the Makefile FILES list: stock -> patched
OBJ_SWAPS = [("coarea.o", "coarea_fast.o"),
             ("hcapen.o", "hcapen_fast.o"),
             ("holcal.o", "holcal_par.o"),
             ("holeen.o", "holeen_par.o"),
             ("sphqpu.o", "sphqpu_par.o"),
             ("h2dmap.o", "h2dmap_par.o"),
             # tsatr_fast reads VMDHole's packed binary coordinate record when the
             # coord file ends .vhb, and is byte-for-byte the stock reader otherwise.
             ("tsatr.o", "tsatr_fast.o"),
             # concal_par's only behavioural change from stock: the "initial
             # point probe radius ... less than probe radius" notice can be
             # RECORDED instead of printed. Unguarded it fires from inside
             # holcal_par's OpenMP prepass, so its position in the .out
             # depended on thread scheduling (measured differing between 1 and
             # 8 threads). The prepass sets CN_REC (a threadprivate common),
             # CONCAL then stores the radius rather than writing, and
             # holcal_par replays the recorded notices in ascending plane
             # order. Gating on SHORTO instead does NOT work: the value comes
             # from a HOLEEN call inside CONCAL, and pass 2 does not re-run
             # CONCAL for prepass-covered planes, so 6 of 13 notices were lost.
             ("concal.o", "concal_par.o")]
# Files compiled with -fopenmp (per-file, never globally - a global -fopenmp puts
# hole.f's large arrays on the stack and segfaults every mode).
# sphqpu_par MUST be here: its parallel dot-culling loop is guarded by !$OMP
# sentinels, so without -fopenmp they are just comments and sph_process silently
# builds SERIAL - correct output, but none of the speedup.
# h2dmap_par: the 2DMAPS wall-distance loop, which the unrolled map waits on.
OMP_FILES = ["holcal_par", "holeen_par", "concal_par", "coarea_fast", "sphqpu_par",
             "h2dmap_par"]


def backup_once(path):
    b = path + ".vmdhole_orig"
    if not os.path.exists(b):
        shutil.copy2(path, b)


def patch_machine_dep(src):
    path = os.path.join(src, "machine_dep.g77")
    if not os.path.isfile(path):
        sys.exit("ERROR: %s not found - is this a hole2/src tree?" % path)
    txt = open(path).read()
    if "RNGCNT" in txt:
        print("  machine_dep.g77 : already patched (skip)")
        return
    orig = txt
    # 1) declare the counter right after the CSEED common block (inside dRAND)
    txt, n1 = re.subn(
        r"(\n\s*COMMON\s*/CSEED/\s*FSEED\n)",
        r"\1"
        "C VMDHole parallel-CONNOLLY: inert RNG draw counter (see holcal_par.f)\n"
        "      INTEGER           NDRAW\n"
        "      COMMON /RNGCNT/   NDRAW\n",
        txt, count=1)
    # 2) zero it on the very first dRAND call (right after LFIRST = .TRUE.)
    txt, n2 = re.subn(
        r"(\n\s*LFIRST\s*=\s*\.TRUE\.\n)",
        r"\1        NDRAW = 0\n",
        txt, count=1)
    # 3) increment it on every draw (the RAND(ISEED) just before dRAND's RETURN)
    txt, n3 = re.subn(
        r"(\n\s*ISEED\s*=\s*0\n\s*RANDOM\s*=\s*RAND\(\s*ISEED\s*\)\n)",
        r"\1      NDRAW = NDRAW + 1\n",
        txt, count=1)
    if not (n1 and n2 and n3):
        sys.exit("ERROR: machine_dep.g77 does not match the expected stock dRAND "
                 "(anchors %d/%d/%d found). Refusing to patch a non-stock file."
                 % (n1, n2, n3))
    backup_once(path)
    open(path, "w").write(txt)
    print("  machine_dep.g77 : NDRAW counter added")


def patch_makefile(src):
    path = os.path.join(src, "Makefile")
    if not os.path.isfile(path):
        sys.exit("ERROR: %s not found" % path)
    txt = open(path).read()
    # PER-ITEM idempotency, NOT a single global marker. A global
    # "already patched -> return" guard silently freezes a tree at whatever the
    # patch set was when it was FIRST patched: a later-added swap (this is exactly
    # how sphqpu_par was missed - the tree kept building a SERIAL sph_process, with
    # correct output and none of the speedup, and re-running this script did
    # nothing) can never be applied. Each swap and each rule is now checked and
    # added on its own, so re-running always converges on the full patch set.
    changed = False
    # object-file swaps (anchored so re-running can't double-swap)
    for stock, fast in OBJ_SWAPS:
        if re.search(r"(?m)^" + re.escape(fast) + r"(\s*\\?)\s*$", txt):
            continue                      # this one is already swapped
        txt, n = re.subn(r"(?m)^" + re.escape(stock) + r"(\s*\\?)\s*$",
                         fast + r"\1", txt)
        if n == 0:
            sys.exit("ERROR: Makefile FILES list has neither '%s' nor '%s' - not a "
                     "stock hole2 Makefile." % (stock, fast))
        changed = True
    # define OMPFLAGS + link with -fopenmp, right after the FFLAGS line
    if "OMPFLAGS := -fopenmp" in txt:
        n = 1                              # already present from an earlier run
    else:
      txt, n = re.subn(
        r"(?m)^(FFLAGS\s*\+=.*\n)",
        r"\1"
        "# VMDHole parallel CONNOLLY: -fopenmp applied ONLY to the 4 parallel-\n"
        "# region files below (per-file rules), never globally - global -fopenmp\n"
        "# implies -frecursive, pushing hole.f's ~24MB arrays onto the stack and\n"
        "# seg-faulting every mode.\n"
        "OMPFLAGS := -fopenmp\n"
        "LFLAGS += -fopenmp\n",
        txt, count=1)
      changed = changed or n > 0
    if n == 0:
        sys.exit("ERROR: Makefile has no 'FFLAGS +=' line.")
    # per-file -fopenmp compile rules, appended at end
    missing = [f for f in OMP_FILES
               if not re.search(r"(?m)^" + re.escape(f) + r"\.o:\s", txt)]
    if missing:
        rules = ["", "# VMDHole parallel CONNOLLY: per-thread automatic arrays +",
                 "# OpenMP directives + THREADPRIVATE for the parallel-region files."]
        for f in missing:
            rules += ["%s.o: %s.f" % (f, f),
                      "\t$(FC) $(FFLAGS) $(OMPFLAGS) -c %s.f" % f,
                      "\t$(AR) -rcv $(LIB_NAME) %s.o" % f, ""]
        txt = txt.rstrip("\n") + "\n\n" + "\n".join(rules) + "\n"
        changed = True
    if not changed:
        print("  Makefile        : already fully patched (skip)")
        return
    backup_once(path)
    open(path, "w").write(txt)
    print("  Makefile        : swaps=%s  omp-rules-added=%s"
          % ([f for _, f in OBJ_SWAPS], missing or "none"))


def copy_new_files(src, here):
    for f in NEW_FILES:
        s = os.path.join(here, f)
        if not os.path.isfile(s):
            sys.exit("ERROR: patch file %s missing next to this script." % f)
        shutil.copy2(s, os.path.join(src, f))
        print("  %-16s: copied" % f)


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: apply_patches.py <hole2/src dir>")
    src = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    if not os.path.isfile(os.path.join(src, "Makefile")):
        sys.exit("ERROR: '%s' is not a hole2/src tree (no Makefile)." % src)
    print(">> Applying VMDHole parallel-CONNOLLY patches to %s" % src)
    check_base_files(src)
    copy_new_files(src, here)
    patch_machine_dep(src)
    patch_makefile(src)
    print(">> Patches applied (originals backed up to *.vmdhole_orig).")


if __name__ == "__main__":
    main()
