#!/bin/sh
# An engine path set BEFORE init_executables must survive it.
#
# THE DEFECT THIS GUARDS:
#   docs/scripting.md tells batch authors to point the plugin at their engines
#   with `set ::VMDHole::state(hole_exec) /path/to/hole` and then call
#   init_executables. init_executables calls load_config FIRST, and load_config
#   assigns every persisted key unconditionally - hole_exec among them. So a
#   ~/.vmdhole_config holding an empty or stale hole_exec (one exists for anyone
#   who has opened the GUI once) silently threw the caller's path away.
#   The run then used the embedded Tcl engine, roughly 100x slower, while the
#   script believed it had selected a compiled one.
#
# Pure Tcl: the two procs are lifted out of the shipped file, so the test cannot
# drift from the code it checks. No VMD needed.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/vmdhole/vmdhole.tcl"

echo "preset-exec-paths: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no vmdhole.tcl"; exit 0; }
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM
# a stand-in engine: only its executability is tested here
printf '#!/bin/sh\nexit 0\n' > "$T/hole"; chmod +x "$T/hole"
# a config that persists an EMPTY hole_exec - the shape that triggers the bug
printf 'hole_exec = \nsph_process_exec = \nsos_triangle_exec = \n' > "$T/cfg"

# A Windows-shaped engine dir: bare "hole"/"sph_process" never exist there,
# only the .exe-suffixed builds native/build.sh actually produces.
mkdir -p "$T/winenv"
printf '#!/bin/sh\nexit 0\n' > "$T/winenv/hole.exe";        chmod +x "$T/winenv/hole.exe"
printf '#!/bin/sh\nexit 0\n' > "$T/winenv/sph_process.exe"; chmod +x "$T/winenv/sph_process.exe"
printf '#!/bin/sh\nexit 0\n' > "$T/winenv/sos_triangle.exe"; chmod +x "$T/winenv/sos_triangle.exe"
printf '#!/bin/sh\nexit 0\n' > "$T/winenv/mole_tunnel_engine.exe"; chmod +x "$T/winenv/mole_tunnel_engine.exe"
# hydro_project must answer --hole-features with "hydroproject" - the proc
# probes it to rule out an unrelated same-named executable.
printf '#!/bin/sh\necho hydroproject\n' > "$T/winenv/hydro_project.exe"; chmod +x "$T/winenv/hydro_project.exe"
printf 'hole_exec = \nsph_process_exec = \nsos_triangle_exec = \n' > "$T/wincfg"

cat > "$T/drv.tcl" <<'TCLEOF'
namespace eval ::VMDHole {}
set SRC [lindex $argv 0]; set CFG [lindex $argv 1]; set ENG [lindex $argv 2]
set WINCFG [lindex $argv 3]; set WINENV [lindex $argv 4]
set fh [open $SRC r]; set src [read $fh]; close $fh

proc lift {src name} {
    set pat "proc ::VMDHole::$name "
    set i [string first $pat $src]
    if {$i < 0} { return "" }
    set buf ""
    foreach ln [split [string range $src $i end] "\n"] {
        append buf $ln "\n"
        if {[info complete $buf]} { return $buf }
    }
    return ""
}
foreach p {load_config _config_skip_keys init_executables find_hole_exe _find_exe \
           _mole_engine_path _hydro_project_path save_config _note} {
    set b [lift $src $p]
    if {$b eq ""} { puts "FATAL: could not lift ::VMDHole::$p"; exit 3 }
    namespace eval ::VMDHole $b
}
# stubs for the bits that need VMD / would write to $HOME
proc vmdcon {args} {}
proc ::VMDHole::save_config {} {}
proc ::VMDHole::_sweep_stale_tmpdirs {} {}
proc ::VMDHole::_note {args} {}

namespace eval ::VMDHole {
    variable config_file
    variable state
    variable _swept 0
}
set ::VMDHole::config_file $CFG
foreach k {hole_exec sph_process_exec sos_triangle_exec mole_engine_exec radius_file} {
    set ::VMDHole::state($k) ""
}

# The documented recipe: set the path, then init.
set ::VMDHole::state(hole_exec) $ENG
::VMDHole::init_executables

if {$::VMDHole::state(hole_exec) eq $ENG} {
    puts "OK preset hole_exec survived init_executables"
} else {
    puts "BAD preset hole_exec was discarded (now '$::VMDHole::state(hole_exec)')"
}

# Guard the guard: the config must really carry an empty hole_exec, or this
# fixture would pass even on the buggy code.
set fh [open $CFG r]; set cfgtext [read $fh]; close $fh
if {[regexp {hole_exec\s*=\s*$} [string trim $cfgtext "\n"]] || [string match "*hole_exec = \n*" $cfgtext]} {
    puts "OK fixture config really persists an empty hole_exec"
} else {
    puts "BAD fixture config no longer triggers the overwrite"
}

# --- .exe discovery: bare names never exist on Windows, only *.exe does ----
foreach k {hole_exec sph_process_exec sos_triangle_exec mole_engine_exec radius_file} {
    set ::VMDHole::state($k) ""
}
set ::VMDHole::config_file $WINCFG
set ::env(VMDHOLE_HOLE_EXE_DIR) $WINENV
::VMDHole::init_executables

set wantHole [file join $WINENV hole.exe]
set wantSph  [file join $WINENV sph_process.exe]
if {$::VMDHole::state(hole_exec) eq $wantHole} {
    puts "OK .exe discovery found hole.exe via VMDHOLE_HOLE_EXE_DIR"
} else {
    puts "BAD .exe discovery: hole_exec = '$::VMDHole::state(hole_exec)', wanted '$wantHole'"
}
if {$::VMDHole::state(sph_process_exec) eq $wantSph} {
    puts "OK .exe discovery backfilled the sph_process.exe sibling"
} else {
    puts "BAD .exe sibling backfill: sph_process_exec = '$::VMDHole::state(sph_process_exec)', wanted '$wantSph'"
}

# --- the two discovery walks OUTSIDE init_executables -----------------------
# _mole_engine_path's fallback beside sos_triangle: clear the state hit that
# init_executables just backfilled, so the sibling WALK itself is what runs.
set ::VMDHole::state(mole_engine_exec) ""
set ::VMDHole::state(sos_triangle_exec) [file join $WINENV sos_triangle.exe]
set wantMole [file join $WINENV mole_tunnel_engine.exe]
set gotMole [::VMDHole::_mole_engine_path]
if {$gotMole eq $wantMole} {
    puts "OK _mole_engine_path's sibling walk found mole_tunnel_engine.exe"
} else {
    puts "BAD _mole_engine_path sibling walk: got '$gotMole', wanted '$wantMole'"
}

# _hydro_project_path: no init_executables backfill exists for it, so this
# walk is its only non-PATH discovery. It also execs the found binary with
# --hole-features and requires "hydroproject" in the reply (the fixture
# answers that). Unset the memo first or a previous "" would be returned.
catch {unset ::VMDHole::_hydro_project_exe}
set wantHydro [file join $WINENV hydro_project.exe]
set gotHydro [::VMDHole::_hydro_project_path]
if {$gotHydro eq $wantHydro} {
    puts "OK _hydro_project_path found hydro_project.exe beside sos_triangle"
} else {
    puts "BAD _hydro_project_path: got '$gotHydro', wanted '$wantHydro'"
}
TCLEOF

out=$(timeout 60 tclsh "$T/drv.tcl" "$SRC" "$T/cfg" "$T/hole" "$T/wincfg" "$T/winenv" 2>&1); rc=$?
printf '%s\n' "$out" | sed 's/^/    /'
if [ "$rc" -eq 3 ] || printf '%s\n' "$out" | grep -q '^FATAL'; then
    echo "  FAIL  the check could not run (see FATAL above)"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
nok=$(printf '%s\n' "$out" | grep -c '^OK ')
nbad=$(printf '%s\n' "$out" | grep -c '^BAD ')
if [ "$nbad" -eq 0 ]; then
    echo "  PASS  a caller-set engine path survives init_executables, and .exe discovery works ($nok checks)"
    echo "  -> 1 passed, 0 failed"
else
    echo "  FAIL  $nbad check(s) failed"
    echo "  -> 0 passed, 1 failed"; exit 1
fi
