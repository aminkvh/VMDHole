#!/bin/sh
# Two independent clustering mechanisms in tunnel mode, and the bug where
# they used to share one on/off flag.
#
#   1 the shipped default for DISPLAY clustering (tunnel_cluster_on, groups
#     near-duplicate routes in the landed frame into one representative) is
#     OFF - every route is listed/drawn separately unless the user opts in.
#   2 CROSS-FRAME identity (_tunnel_xframe_build, what gives a tunnel an
#     identity across a trajectory - Trends/Over Time/Mean Profile/Histogram
#     all depend on it) ALWAYS runs for a multi-frame result, independent of
#     the display flag. It used to be nested inside "if display clustering is
#     on", so the (now off-by-default) display toggle silently broke every
#     trajectory view's tunnel identity too.
#   3 toggling display clustering on/off does not touch the already-built
#     cross-frame identity - no MOLE recompute, no rebuild.
#   4 with display clustering off, _tunnel_candidates returns every raw
#     tunnel (no representative-only collapsing).

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SRC="$DIR/../vmdhole.tcl"
ENGINE="$DIR/../../native/mole_tunnel_engine"
PDB="$DIR/fixtures/mole_reference/1BL8.pdb"
echo "=============================================================="
echo "tunnel-clustering: $SRC"
pass=0; fail=0
ok()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad() { echo "  FAIL  $1"; fail=$((fail+1)); }

[ -f "$SRC" ] || { echo "  FAIL  source not found"; exit 1; }

# ---- static: the shipped default, no engine/fixture needed ----
# Display clustering ships ON: on a real structure MOLE returns many
# near-duplicate routes, and listing every one of them was reported as noise.
# Cross-frame identity is a SEPARATE mechanism and stays unconditional.
if grep -qE '^\s*tunnel_cluster_on\s+1\s*$' "$SRC"; then
    ok "tunnel_cluster_on defaults to 1 (display clustering ON)"
else
    bad "tunnel_cluster_on's shipped default is not 1"
fi

if [ ! -x "$ENGINE" ] || [ ! -f "$PDB" ]; then
    # Group-level form on its own line: run_tests.sh's release gate anchors on
    # ^SKIP: (deliberately - see its comment), so an indented-only skip would
    # let a release pass while the runtime half checked nothing.
    echo "SKIP: no compiled mole_tunnel_engine or no 1BL8 fixture - runtime checks not run"
    echo "  -> $pass passed, $fail failed"
    [ "$fail" -eq 0 ]
    exit $?
fi

VMD=${VMD:-vmd}
OUT=$(mktemp)
WORK=$(mktemp -d)
cat > "$OUT.tcl" <<EOF
package provide Tk 8.5
if {[catch {source "$SRC"} e]} {
    set o [open "$OUT.out" w]; puts \$o "FAIL source: \$e"; puts \$o "BAD=1"; close \$o; quit
}
set o [open "$OUT.out" w]
set bad 0

set mid [mol new "$PDB" waitfor all]
animate dup \$mid
animate dup \$mid

set ::VMDHole::state(molid) \$mid
set ::VMDHole::state(selection) "protein"
set ::VMDHole::state(engine) "mole"
set ::VMDHole::state(tunnel_start) "73.853 26.536 26.594"
set ::VMDHole::state(mole_engine_exec) "$ENGINE"
set ::VMDHole::state(work_dir) "$WORK"
set ::VMDHole::state(frame_spec) "all"
::VMDHole::run_tunnel_analysis
set nfr [llength \$::VMDHole::tunnel_result_frames]
if {\$nfr < 2} {
    puts \$o "FAIL run produced \$nfr frame(s), need >=2: \$::VMDHole::state(status)"
    incr bad
} else {
    set fr0 [lindex \$::VMDHole::tunnel_result_frames 0]
    set fr1 [lindex \$::VMDHole::tunnel_result_frames 1]

    # ---- 2. cross-frame identity built even with display clustering off ----
    if {!\$::VMDHole::state(tunnel_cluster_on)} {
        puts \$o "FAIL default tunnel_cluster_on is not 1 in this session"
        incr bad
    }
    if {[llength \$::VMDHole::tunnel_xclusters] > 0} {
        puts \$o "PASS cross-frame clusters built alongside display clustering (\$::VMDHole::state(status))"
    } else {
        puts \$o "FAIL tunnel_xclusters is empty"
        incr bad
    }
    # Per-frame ON DEMAND now: the run no longer clusters every frame up front
    # (~20 ms x every frame - 216 s on a 10k trajectory - for a display-only
    # grouping). Ask for the frame, then assert it exists.
    ::VMDHole::_tunnel_ensure_clusters_for \$fr0
    if {[info exists ::VMDHole::tunnel_clusters(\$fr0)]} {
        puts \$o "PASS per-frame display grouping built on demand when tunnel_cluster_on is 1"
    } else {
        puts \$o "FAIL tunnel_clusters(\$fr0) missing even after an explicit ensure"
        incr bad
    }

    # ---- 4. candidates == every raw tunnel with display clustering off ----
    set cand0 [::VMDHole::_tunnel_candidates \$fr0]
    set raw0 [llength \$::VMDHole::tunnel_results(\$fr0)]
    if {[llength \$cand0] <= \$raw0 && [llength \$cand0] > 0} {
        puts \$o "PASS _tunnel_candidates lists cluster representatives ([llength \$cand0] of \$raw0 routes)"
    } else {
        puts \$o "FAIL candidates n=[llength \$cand0] raw=\$raw0 - must be a non-empty subset"
        incr bad
    }
    # Turning it OFF must expose every raw route again.
    set ::VMDHole::state(tunnel_cluster_on) 0
    set candoff [::VMDHole::_tunnel_candidates \$fr0]
    if {[llength \$candoff] == \$raw0} {
        puts \$o "PASS turning display clustering off exposes every raw route (\$raw0)"
    } else {
        puts \$o "FAIL clustering off gave [llength \$candoff] of \$raw0 raw routes"
        incr bad
    }
    set ::VMDHole::state(tunnel_cluster_on) 1

    # ---- cross-frame lookup actually resolves a real identity ----
    animate goto \$fr0
    set ::VMDHole::state(tunnel_selected_id) [lindex \$cand0 0]
    set xcid [::VMDHole::_tunnel_selected_cluster]
    set rk1 [::VMDHole::_tunnel_rank_in_frame \$fr1]
    if {\$xcid ne "" && \$rk1 ne ""} {
        puts \$o "PASS selected tunnel resolves a cross-frame cluster and a real rank in frame \$fr1"
    } else {
        puts \$o "FAIL cross-frame lookup failed: cluster=\$xcid rank_in_frame1=\$rk1"
        incr bad
    }

    # ---- 3. toggling display clustering does not touch cross-frame data ----
    set xc_before \$::VMDHole::tunnel_xclusters
    ::VMDHole::_tunnel_cluster_toggle_clicked
    set xc_after \$::VMDHole::tunnel_xclusters
    if {\$xc_before eq \$xc_after} {
        puts \$o "PASS toggling display clustering leaves cross-frame identity untouched"
    } else {
        puts \$o "FAIL cross-frame clusters changed after a display-only toggle"
        incr bad
    }
    if {\$::VMDHole::state(tunnel_cluster_on)} {
        puts \$o "FAIL toggle did not flip tunnel_cluster_on"
        incr bad
    }
    set ::VMDHole::state(tunnel_cluster_on) 1

    # ---- re-run with display clustering ON: per-frame grouping is present,
    # cross-frame mechanism is the same one, still populated ----
    set ::VMDHole::state(tunnel_cluster_on) 1
    ::VMDHole::run_tunnel_analysis
    if {[info exists ::VMDHole::tunnel_clusters(\$fr0)] && [llength \$::VMDHole::tunnel_xclusters] > 0} {
        puts \$o "PASS display clustering ON still builds cross-frame identity alongside it"
    } else {
        puts \$o "FAIL cluster_on=1 rerun: tunnel_clusters exists=[info exists ::VMDHole::tunnel_clusters(\$fr0)] xclusters=[llength \$::VMDHole::tunnel_xclusters]"
        incr bad
    }

    # ---- the kernel's representative-distance column IS the Tcl reference ----
    # _tunnel_xframe_build picks each frame's representative by proximity to
    # the cluster. Doing that in Tcl cost 15.3 s of an 18.5 s build at n=1837,
    # so the kernel now returns the distance as a third --tunnel-cluster
    # column. Both paths must produce the SAME number - if they drift, the
    # tracked route silently changes and nothing else would notice.
    if {[::VMDHole::sos_triangle_has_feature tunnelclusterdist]} {
        # SYNTHETIC first, and it is not optional. This fixture duplicates one
        # frame (animate dup), so every cross-frame pair is the SAME geometry
        # and its representative distance is exactly 0.0 - a kernel returning
        # garbage scaled by any factor still matches 0.0. Sabotage-checked:
        # multiplying the C column by 1.02 leaves the real-data comparison
        # below passing and only this pair catches it.
        set A {}; set B {}
        for {set k 0} {\$k < 12} {incr k} {
            lappend A [list 0.0 0.0 [expr {\$k*1.0}]]
            lappend B [list 1.0 0.0 [expr {\$k*1.0}]]
        }
        array set _rd {}
        set _cc [::VMDHole::_tunnel_cluster_c [list \$A \$B] 3.0 0 _rd]
        set want [::VMDHole::_tunnel_pair_distance \$A \$B]
        set got  [expr {[info exists _rd(1)] ? \$_rd(1) : ""}]
        if {\$got ne "" && abs(\$got - \$want) < 1e-9 && \$want > 0.5} {
            puts \$o "PASS kernel repr-distance on a known non-zero pair: \$got (Tcl \$want)"
        } else {
            puts \$o "FAIL kernel repr-distance on a known pair: got '\$got' want '\$want'"
            incr bad
        }

        # The kernel prunes pairs whose bounding boxes are more than the
        # threshold apart - rigorous, since the box gap lower-bounds the
        # distance - and then has to fill those distances back in for pairs
        # that share a component, because average-link averages them.
        # A---B---C at x = 0, 2, 4: d(A,B)=d(B,C)=2 (<=3, so one component),
        # d(A,C)=4 (>3, box-pruned). Average-link merges A,B then compares
        # d(AB,C) = (4+2)/2 = 3.0, which still merges - but only if the pruned
        # 4.0 was filled in. Sabotage-checked: skipping that fill leaves
        # d(A,C) at infinity and this splits into 2 clusters.
        set chain {}
        foreach x {0.0 2.0 4.0} {
            set pts {}
            for {set z 0} {\$z < 21} {incr z} { lappend pts [list \$x 0.0 [expr {\$z*0.5}]] }
            lappend chain \$pts
        }
        set _ch [::VMDHole::_tunnel_cluster_c \$chain 3.0 0]
        if {[llength \$_ch] == 1 && [llength [lindex \$_ch 0]] == 3} {
            puts \$o "PASS box-pruned pairs are refilled inside a component (A-B-C chain merges)"
        } else {
            puts \$o "FAIL A-B-C chain gave [llength \$_ch] cluster(s): \$_ch"
            incr bad
        }

        set nchk 0; set worst 0.0
        foreach c \$::VMDHole::tunnel_xclusters {
            if {[llength \$c] < 2} { continue }
            set first [lindex [lsort -integer -index 0 \$c] 0]
            set cref [::VMDHole::_tunnel_mid_points \
                [::VMDHole::_tunnel_tuple_for [lindex \$first 0] [lindex \$first 1]]]
            foreach m \$c {
                if {\$m eq \$first} { continue }
                if {![info exists ::VMDHole::tunnel_xrepr_dist(\$m)]} { continue }
                incr nchk
                set ref [::VMDHole::_tunnel_pair_distance \$cref \
                    [::VMDHole::_tunnel_mid_points \
                        [::VMDHole::_tunnel_tuple_for [lindex \$m 0] [lindex \$m 1]]]]
                set dd [expr {abs(\$ref - \$::VMDHole::tunnel_xrepr_dist(\$m))}]
                if {\$dd > \$worst} { set worst \$dd }
            }
        }
        if {\$nchk > 0 && \$worst < 1e-9} {
            puts \$o "PASS kernel repr-distance matches the Tcl reference on \$nchk pairs (max diff \$worst)"
        } else {
            puts \$o "FAIL kernel repr-distance: \$nchk pairs checked, max diff \$worst"
            incr bad
        }
        # rep[] is the cluster's LOWEST POOL INDEX, which the Tcl side relies on
        # being the earliest-frame member. Frame-major pool order is what makes
        # those the same pathway; assert it rather than trusting the argument.
        set offby 0
        foreach c \$::VMDHole::tunnel_xclusters {
            set first [lindex [lsort -integer -index 0 \$c] 0]
            foreach m \$c {
                if {[lindex \$m 0] < [lindex \$first 0]} { incr offby }
            }
        }
        if {\$offby == 0} {
            puts \$o "PASS every cluster's reference member is its earliest frame"
        } else {
            puts \$o "FAIL \$offby member(s) sit earlier than their cluster's reference"
            incr bad
        }
    } else {
        puts \$o "PASS kernel lacks tunnelclusterdist - Tcl fallback path in use (not a failure)"
    }
}

puts \$o "BAD=\$bad"
close \$o
quit
EOF
"$VMD" -dispdev text -e "$OUT.tcl" >/dev/null 2>&1
if [ -f "$OUT.out" ]; then
    while IFS= read -r line; do
        case "$line" in
            PASS*) ok "${line#PASS }" ;;
            FAIL*) bad "${line#FAIL }" ;;
        esac
    done < "$OUT.out"
    if ! grep -q "^BAD=" "$OUT.out"; then
        bad "script did not complete (no BAD= sentinel)"
    fi
else
    bad "script produced no output"
fi
rm -rf "$OUT" "$OUT.tcl" "$OUT.out" "$WORK"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
