#!/bin/sh
# Tunnel mode must never write into, read from, or overwrite HOLE's results.
#
# The invariants below are the ones a future change can break WITHOUT anything
# failing loudly: both modes write frame-indexed folders, and both feed the same
# surface pipeline, so a collision shows up as "my HOLE results turned into
# tunnels" long after the commit that caused it. Each check here maps to one
# decision in feat(tunnel) wiring:
#
#   1 output roots differ, and an assigned work_dir puts tunnels in a subfolder
#     rather than a sibling, so "all my files in the folder I assigned" holds.
#   2 tunnel frame folders are tunnel_%05d, so collect_frame_dirs - which drives
#     HOLE's Import - cannot pick them up.
#   3 run_tunnel_analysis does not touch last_root_dir / import_dir. Those aim
#     HOLE's Import; pointing them at tunnel output is the exact confusion the
#     separation exists to prevent. Static, because the runtime path needs a
#     structure with a real cavity.
#   4 Pore Profile/Trends are now deliberately SHARED (a selected MOLE tunnel
#     has its own Radius-vs-distance profile), so this checks that the
#     tunnel-specific draw procs never reference HOLE's own result storage,
#     and that the shared draw procs branch on analysis_mode before touching
#     it - data independence, not tab identity, is what actually matters.

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SRC="$DIR/../vmdhole.tcl"
echo "=============================================================="
echo "tunnel-separation: $SRC"
pass=0; fail=0
ok()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad() { echo "  FAIL  $1"; fail=$((fail+1)); }

[ -f "$SRC" ] || { echo "  FAIL  source not found"; exit 1; }

# ---- 3. static: the import-state variables are never written by tunnel mode ----
body=$(awk '/^proc ::VMDHole::run_tunnel_analysis /,/^}/' "$SRC")
if printf '%s' "$body" | grep -qE 'set +state\((last_root_dir|import_dir)\)'; then
    bad "run_tunnel_analysis assigns last_root_dir/import_dir (drives HOLE's Import)"
else
    ok "run_tunnel_analysis leaves last_root_dir / import_dir alone"
fi
if printf '%s' "$body" | grep -qE 'format +"frame_%0'; then
    bad "run_tunnel_analysis uses frame_%05d - collect_frame_dirs would match it"
else
    ok "run_tunnel_analysis does not use HOLE's frame_%05d naming"
fi

# ---- 1, 2, 4: runtime ----
VMD=${VMD:-vmd}
# Static checks above ran under plain tclsh; the rest drives a real VMD.
# SKIP (not FAIL) without one - see test_tcl_pitfalls.sh for the convention.
if ! command -v "$VMD" >/dev/null 2>&1; then
    echo "  SKIP  runtime separation checks (no vmd on PATH)"
    echo "  -> $pass passed, $fail failed"
    [ "$fail" -eq 0 ]; exit $?
fi
OUT=$(mktemp)
WORK=$(mktemp -d)
cat > "$OUT.tcl" <<'EOF'
package provide Tk 8.5
if {[catch {source $::env(SEP_SRC)} e]} {
    set o [open $::env(SEP_OUT) w]; puts $o "FAIL source: $e"; puts $o "BAD=1"; close $o; quit
}
set o [open $::env(SEP_OUT) w]
set bad 0
set work $::env(SEP_WORK)

# A molecule is needed only for get_molecule_basename; molid 0 with no mol is
# fine because work_dir is assigned, which short-circuits before that lookup.
set ::VMDHole::state(work_dir) $work
set hole   [lindex [::VMDHole::resolve_output_root 0]        0]
set tunnel [lindex [::VMDHole::resolve_output_root 0 tunnel] 0]
if {$hole eq $tunnel} { puts $o "FAIL roots identical: $hole"; incr bad }
if {$tunnel ne [file join $hole tunnels]} {
    puts $o "FAIL tunnel root not <work_dir>/tunnels: $tunnel"; incr bad
}

# collect_frame_dirs must see HOLE's folder and ignore the tunnel one.
file mkdir [file join $work frame_00000]
file mkdir [file join $work tunnel_00000]
set found [::VMDHole::collect_frame_dirs $work]
if {[llength $found] != 1} { puts $o "FAIL collect_frame_dirs found [llength $found]: $found"; incr bad }
foreach f $found {
    if {[string match "*tunnel_*" $f]} { puts $o "FAIL collect_frame_dirs matched a tunnel dir: $f"; incr bad }
}
set tfound [::VMDHole::collect_frame_dirs [file join $work tunnels]]
if {[llength $tfound] != 0} { puts $o "FAIL tunnel root exposes frame dirs: $tfound"; incr bad }

# Pore Profile and Trends are now deliberately SHARED between modes (a single
# selected MOLE tunnel has its own well-defined Radius-vs-distance profile,
# the same shape HOLE's own Pore Profile plots) - draw_profile_plot and
# draw_minr_tab branch on analysis_mode into draw_tunnel_profile_plot/
# draw_tunnel_trends_plot BEFORE touching anything HOLE-specific. What must
# still hold, and is checked here in its place, is that the tunnel-mode draw
# procs never reference HOLE's own result storage - that is the actual
# separation invariant now, not tab identity.
set tt [::VMDHole::_mode_tab_set tunnel]
if {[llength $tt] == 0} { puts $o "FAIL tunnel mode shows no tabs"; incr bad }
foreach p {draw_tunnel_profile_plot draw_tunnel_trends_plot} {
    set body [info body ::VMDHole::$p]
    # results is a DICT (variable results [dict create]), read via
    # [dict get/exists $results ...] - never "results(" array syntax, so the
    # load-bearing check is "does this proc even bring the dict into scope".
    # Without "variable results" a bare $results is just an uninitialized
    # local, not HOLE's data, so this alone is sufficient - not a heuristic.
    foreach forbidden {{variable results} ensure_profile_full} {
        if {[string match "*$forbidden*" $body]} {
            puts $o "FAIL $p references HOLE's own '$forbidden' - the two modes' Pore Profile/Trends must stay data-independent even though they now share tabs"
            incr bad
        }
    }
}
# And the reverse: draw_profile_plot/draw_minr_tab must branch to the tunnel
# proc before doing any HOLE-specific work, not after.
foreach {p want} {draw_profile_plot draw_tunnel_profile_plot draw_minr_tab draw_tunnel_trends_plot} {
    set body [info body ::VMDHole::$p]
    set branch_at [string first "analysis_mode" $body]
    set touch_at  [string first "variable results" $body]
    if {$branch_at < 0} { puts $o "FAIL $p has no analysis_mode branch at all"; incr bad }
    if {$touch_at >= 0 && $branch_at > $touch_at} {
        puts $o "FAIL $p touches HOLE's results before checking analysis_mode"; incr bad
    }
}

# Results storage must be distinct variables, not aliases.
set ::VMDHole::result_frames {7}
set ::VMDHole::tunnel_result_frames {}
if {[llength $::VMDHole::result_frames] != 1 || [llength $::VMDHole::tunnel_result_frames] != 0} {
    puts $o "FAIL result_frames and tunnel_result_frames are not independent"; incr bad
}

puts $o "BAD=$bad"
close $o
quit
EOF
SEP_OUT="$OUT.out" SEP_SRC="$SRC" SEP_WORK="$WORK" $VMD -dispdev text -e "$OUT.tcl" >/dev/null 2>&1
if [ -f "$OUT.out" ] && grep -q "^BAD=0$" "$OUT.out"; then
    ok "roots differ, work_dir/tunnels honoured, Import blind to tunnel dirs, shared tabs stay data-independent"
else
    bad "a separation invariant broke:"
    [ -f "$OUT.out" ] && sed 's/^/          /' "$OUT.out"
fi
rm -rf "$OUT" "$OUT.tcl" "$OUT.out" "$WORK"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
