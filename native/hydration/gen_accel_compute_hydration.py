#!/usr/bin/env python3
"""Generate _accel_compute_hydration: an accelerated variant of
compute_hydration that calls native/hydro_project for the measured
Phase A (COG reduction + axial projection + envelope-radius test) and Phase C
(KDE/hard binning) work, while extracting every OTHER line of the proc
VERBATIM (byte-identical) from vmdhole.tcl via the same anchor-verified,
offset-based introspection as gen_instrumented.py -- read-only, never edits
vmdhole.tcl. Only the water-COG/projection/binning sub-block is hand-written
new Tcl; everything else (bulk density, axis resolution, envelope/bbox
building, Welford aggregation, energy assembly, dict output, file I/O) is the
UNCHANGED original text, so there is no hand-transcription risk in the parts
that are supposed to be identical.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, "..", "..", "vmdhole", "vmdhole.tcl"))
OUT = os.path.join(HERE, "accel_compute_hydration.tcl")

with open(SRC) as f:
    lines = f.readlines()

start_idx = None
for i, l in enumerate(lines):
    if l.startswith("proc ::VMDHole::compute_hydration {} {"):
        start_idx = i
        break
if start_idx is None:
    sys.exit("FATAL: could not find compute_hydration proc start")

# The proc's END is NOT a fixed offset: this file is edited concurrently by
# other agents/sessions, and the tail (Phase D onward, extracted verbatim)
# has already been observed to shift by a few lines between two runs in the
# same session. Find it dynamically: the next top-level "proc ::VMDHole::"
# definition marks the end of this one; walk back from there past any blank
# lines to the lone "}" that closes compute_hydration.
import re as _re
next_proc_idx = None
for i in range(start_idx + 1, len(lines)):
    if _re.match(r'^proc ::VMDHole::\w', lines[i]):
        next_proc_idx = i
        break
if next_proc_idx is None:
    sys.exit("FATAL: could not find the next proc after compute_hydration")
end_idx = next_proc_idx - 1
while end_idx > start_idx and lines[end_idx].strip() == '':
    end_idx -= 1
if lines[end_idx].strip() != '}':
    sys.exit("FATAL: expected compute_hydration's closing brace just before %r, got %r"
              % (lines[next_proc_idx], lines[end_idx]))
end_offset = end_idx - start_idx  # relative offset of the closing "}" line

# --- Anchor location: SEQUENTIAL TEXT SEARCH, not fixed offsets. ---
# vmdhole.tcl is edited concurrently (by other agents, and - once wired - by
# a dispatcher prologue added to
# compute_hydration itself, which shifts every line after it). A fixed-
# offset anchor table breaks under ANY line-count change anywhere before the
# anchor, which has already happened mid-session (the proc's tail shrank by
# 3 lines between two runs). Instead, walk forward from start_idx searching
# for each anchor's (verified-unique-at-this-point) substring, each search
# starting where the previous one left off - this tolerates insertions/
# deletions ANYWHERE in the body, not just after the last anchor found.
# Each entry is either a plain substring (matched via `in`) or a callable
# predicate on the raw line - needed for '}' below, which is NOT unique by
# substring: the per-frame foreach body has several inner-block closing
# braces (8-space indent) before its OWN closing brace (4-space indent, same
# as the `foreach` itself). Matching the bare substring would stop at the
# first (wrong, inner) one.
_close_foreach = lambda l: l == '    }\n' or l == '    }'
ANCHOR_TEXTS = [
    'set fi 0',                                                    # Phase A loop start
    'if {[catch {set wa [atomselect $molid $q frame $frame]}',     # atomselect start
    '$wa delete',
    # ORDER MATTERS: the walk is sequential, so these must appear in the same
    # order as the source. `incr nfwater` moved ABOVE the COG comment when the
    # C accelerator (_hydro_qco_c) was added and the Tcl COG block became its
    # `if {$qco eq ""}` fallback - which silently broke this generator until the
    # two entries were swapped to match. Both are POSITION MARKERS only: nothing
    # is extracted between them (the accelerated variant supplies its own
    # WATER_BLOCK for the nfwater line), so their job is purely to keep the walk
    # anchored across the region that gets dropped.
    'if {[llength $wpos] > 0} { incr nfwater }',                   # counted before the COG block
    '# CHAP maps each solvent RESIDUE',                            # COG block start (dropped)
    'dict set frame_ctx $frame $qco',                              # end of original qco-based frame_ctx set
    'dict set frame_range $frame',
    'dict set frame_radii $frame $_frame_radii',
    'if {$fi % 20 == 0',                                           # progress message
    _close_foreach,                                                # close of per-frame foreach (NOT any '}')
    'set fi 0',                                                    # Phase C loop start (2nd occurrence)
    'if {$nfdata == 0}',                                           # Phase C end / nfdata guard
    'set _T $state(water_temp)',                                   # Phase D start
]
anchor_idx = []  # absolute line index of each anchor, in order
search_from = start_idx
for text in ANCHOR_TEXTS:
    found = None
    matcher = text if callable(text) else (lambda l, _t=text: _t in l)
    for i in range(search_from, end_idx + 1):
        if matcher(lines[i]):
            found = i
            break
    if found is None:
        sys.exit("FATAL: anchor text not found (searched from line %d): %r"
                  % (search_from + 1, text))
    anchor_idx.append(found)
    search_from = found + 1
(a_phaseA_start, a_atomselect_start, a_atomselect_end, a_cog_start, a_cog_end,
 a_frame_ctx, a_frame_range, a_frame_radii, a_progress, a_loop_close,
 a_phaseC_start, a_phaseC_end, a_phaseD_start) = anchor_idx

if not any('return 1' in l for l in lines[a_phaseD_start:end_idx + 1]):
    sys.exit("FATAL: no 'return 1' found between Phase D start and the closing brace"
              " - compute_hydration's tail may have changed structurally, not just length")

def extract_idx(lo, hi):
    """Verbatim lines [lo, hi) by ABSOLUTE index."""
    return "".join(lines[lo:hi])

HEADER = extract_idx(start_idx, a_phaseA_start)        # proc decl through just before the per-frame loop
LOOP_HEAD = extract_idx(a_phaseA_start, a_atomselect_start)  # loop start, validate, axis resolve, env/bbox, q
ATOMSELECT = extract_idx(a_atomselect_start, a_atomselect_end + 1)  # atomselect+get wpos/wres+delete (VMD-bound)
FRAME_TAIL = extract_idx(a_frame_range, a_frame_radii + 1)   # frame_range/frame_radii dict set (the AMISE-
                                   # trigger block between frame_ctx's original set and frame_range is
                                   # dropped entirely: bw_auto=0 is enforced by GUARD, so it would be dead
                                   # code, and it references $qco which this accelerated variant never
                                   # computes in Tcl)
LOOP_CLOSE = extract_idx(a_progress, a_phaseC_start)    # progress message, closing brace of the per-frame
                                                          # foreach, blank line, Phase B (dead: frame_qco_
                                                          # need_bw is always empty under the GUARD)
POST_D_ONWARD = extract_idx(a_phaseD_start, end_idx + 1)  # Phase D onward through closing brace, UNCHANGED

HEADER = HEADER.replace(
    'proc ::VMDHole::compute_hydration {} {',
    'proc ::VMDHole::_accel_compute_hydration {} {', 1)

# Guard: this accelerated path only supports a FIXED bandwidth (bw_auto=0).
# bw_auto is computed earlier in HEADER; insert the guard right before the
# per-frame loop begins (same place hp_mark's phaseA_start marker went).
GUARD = (
    '    if {$bw_auto} {\n'
    '        error "hydro_project accelerated path requires a FIXED water_kde_bw'
    ' - set state(water_kde_bw) to a number rather than auto/AMISE"\n'
    '    }\n'
    '    set _hp_batch_dir $::HP_ACCEL_BATCH_DIR\n'
    '    file mkdir $_hp_batch_dir\n'
    '    set _hp_batch_lines {}\n'
    '    set _hp_frame_order {}\n'
)

# Replacement for offsets 289-326 (COG reduction through qco/frame_ctx build):
# write this frame's raw (post-atomselect, pre-COG) water arrays + axis/env
# to a per-frame job file for hydro_project, in the SAME text format
# hydro_project.c parses (parse_and_project() is its definition). nfwater is
# incremented from the RAW wpos length (equivalent to the original's
# post-COG check -- non-empty raw implies non-empty reduced, since COG
# reduction maps at least one atom to at least one residue).
WATER_BLOCK = '''        if {[llength $wpos] > 0} { incr nfwater }
        set _hp_in [file join $_hp_batch_dir [format "f%06d.in" $frame]]
        set _hp_out [file join $_hp_batch_dir [format "f%06d.out" $frame]]
        set _hp_fh [open $_hp_in w]
        puts $_hp_fh [format "AXIS %.17g %.17g %.17g %.17g %.17g %.17g" $mx $my $mz $ux $uy $uz]
        puts $_hp_fh [format "RANGE %.17g %.17g" $cmin $cmax]
        puts $_hp_fh [format "DCAP %.17g" $dcap]
        set _hp_ecos {}; set _hp_ers {}
        foreach _hp_e $env { lassign $_hp_e _hp_co _hp_r; lappend _hp_ecos [format %.17g $_hp_co]; lappend _hp_ers [format %.17g $_hp_r] }
        puts $_hp_fh "ENV [llength $env]"
        puts $_hp_fh [join $_hp_ecos " "]
        puts $_hp_fh [join $_hp_ers " "]
        set _hp_xs {}; set _hp_ys {}; set _hp_zs {}
        foreach _hp_p $wpos { lassign $_hp_p _hp_x _hp_y _hp_z; lappend _hp_xs [format %.17g $_hp_x]; lappend _hp_ys [format %.17g $_hp_y]; lappend _hp_zs [format %.17g $_hp_z] }
        puts $_hp_fh "WATERS [llength $wpos]"
        puts $_hp_fh [join $wres " "]
        puts $_hp_fh [join $_hp_xs " "]
        puts $_hp_fh [join $_hp_ys " "]
        puts $_hp_fh [join $_hp_zs " "]
        puts $_hp_fh [format "BIN %.17g %.17g %d" $dz $bw $use_kde]
        close $_hp_fh
        lappend _hp_batch_lines "$_hp_in\\t$_hp_out"
        lappend _hp_frame_order $frame
        dict set frame_ctx $frame {}
'''

# Replacement for the ENTIRE Phase B+C block (offsets 360-417 in the
# original): run the batch, then reconstruct bincount (global) and
# perframe_raw from hydro_project's output files.
BATCH_CALL = '''
    set frame_bw [dict create]
    set bincount [dict create]
    set _hp_batch_file [file join $_hp_batch_dir "batch.txt"]
    set _hp_global_file [file join $_hp_batch_dir "global.out"]
    set _hp_bf [open $_hp_batch_file w]
    foreach _l $_hp_batch_lines { puts $_hp_bf $_l }
    close $_hp_bf
    if {[llength $_hp_frame_order] > 0} {
        # Support of the KDE sum: every bin the pore's own radius data covers,
        # exactly as the pure-Tcl Phase C computes it. It is global and only
        # becomes known once every frame has been scanned, which is why it goes
        # on the command line rather than in each frame's per-frame BIN line.
        set _hp_kb [lsort -integer [dict keys $binrn]]
        set _hp_klo [expr {[llength $_hp_kb] ? [lindex $_hp_kb 0] : 0}]
        set _hp_khi [expr {[llength $_hp_kb] ? [lindex $_hp_kb end] : -1}]
        set _hp_cmd "[shell_quote $::HP_ACCEL_BIN] --bin --bin-global [shell_quote $_hp_global_file] --kde-range=$_hp_klo,$_hp_khi --batch [shell_quote $_hp_batch_file]"
        if {[catch {exec sh -c $_hp_cmd} _hp_err]} {
            error "hydro_project batch failed: $_hp_err"
        }
        set _hp_gf [open $_hp_global_file r]
        set _hp_gline1 [gets $_hp_gf]
        set _hp_gline2 [gets $_hp_gf]
        close $_hp_gf
        if {[lindex $_hp_gline1 0] ne "BINS"} { error "hydro_project: malformed global output" }
        set _hp_glo [lindex $_hp_gline1 1]; set _hp_ghi [lindex $_hp_gline1 2]
        if {$_hp_ghi >= $_hp_glo} {
            set _hp_gvals $_hp_gline2
            set _hp_bi $_hp_glo
            foreach _hp_v $_hp_gvals { dict set bincount $_hp_bi $_hp_v; incr _hp_bi }
        }
        foreach frame $_hp_frame_order {
            set _hp_out [file join $_hp_batch_dir [format "f%06d.out" $frame]]
            set _hp_fh [open $_hp_out r]
            set _hp_l1 [gets $_hp_fh]
            set _hp_l2 [gets $_hp_fh]
            set _hp_l3 [gets $_hp_fh]
            set _hp_l4 {}
            if {[lindex $_hp_l3 0] eq "BINS"} { set _hp_l4 [gets $_hp_fh] }
            close $_hp_fh
            incr total_w [lindex $_hp_l1 1]
            set frame_bincount [dict create]
            if {[lindex $_hp_l3 0] eq "BINS"} {
                set _hp_lo [lindex $_hp_l3 1]; set _hp_hi [lindex $_hp_l3 2]
                if {$_hp_hi >= $_hp_lo} {
                    set _hp_bi $_hp_lo
                    foreach _hp_v $_hp_l4 { dict set frame_bincount $_hp_bi $_hp_v; incr _hp_bi }
                }
            }
            set _fcmin -1e30; set _fcmax 1e30
            if {[dict exists $frame_range $frame]} { lassign [dict get $frame_range $frame] _fcmin _fcmax }
            set _fr [expr {[dict exists $frame_radii $frame] ? [dict get $frame_radii $frame] : [dict create]}]
            lappend perframe_raw [dict create frame $frame bins $frame_bincount radii $_fr cmin $_fcmin cmax $_fcmax]
        }
    }
'''

body = (HEADER + GUARD + LOOP_HEAD + ATOMSELECT + WATER_BLOCK + FRAME_TAIL
        + LOOP_CLOSE + BATCH_CALL + POST_D_ONWARD)

with open(OUT, "w") as f:
    f.write("# AUTO-GENERATED accelerated variant of compute_hydration (see gen_accel_compute_hydration.py).\n")
    f.write("# Phase A COG+projection and Phase C binning are replaced with calls to\n")
    f.write("# native/hydro_project; every other line is extracted VERBATIM from the\n")
    f.write("# real vmdhole.tcl (never edited). Requires $::HP_ACCEL_BIN and\n")
    f.write("# $::HP_ACCEL_BATCH_DIR to be set before calling.\n")
    f.write(body)
    f.write("\n")

print("OK: wrote", OUT)
