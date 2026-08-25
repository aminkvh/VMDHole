#!/bin/sh
# VMDHole installer. Checks what is required, reports what is missing, then
# copies the plugin into VMD's user plugin directory.
#
#   ./install.sh              install for the current user
#   ./install.sh --check      report requirements and exit, install nothing
#   ./install.sh --dir DIR    install somewhere other than ~/.vmd/plugins
set -e
SRC=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DEST="$HOME/.vmd/plugins"
CHECK_ONLY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --check) CHECK_ONLY=1 ;;
        --dir)   shift; DEST="$1" ;;
        -h|--help) sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done

miss=0
note() { printf '  %-28s %s\n' "$1" "$2"; }

echo "VMDHole requirements"

# VMD itself. Everything else is optional; without VMD there is nothing to load.
VMD=""
for c in vmd vmd2; do command -v "$c" >/dev/null 2>&1 && { VMD=$c; break; }; done
if [ -n "$VMD" ]; then
    note "VMD" "found ($(command -v $VMD))"
else
    note "VMD" "MISSING - required. https://www.ks.uiuc.edu/Research/vmd/"
    miss=1
fi

# HOLE 2. Optional: the plugin ships a pure-Tcl HOLE engine and falls back to it
# automatically, so a missing binary costs speed, not capability.
HOLE=""
for p in "$HOME/hole2/exe/hole" "$(command -v hole 2>/dev/null)"; do
    [ -n "$p" ] && [ -x "$p" ] && { HOLE=$p; break; }
done
if [ -n "$HOLE" ]; then
    note "HOLE 2" "found ($HOLE)"
else
    note "HOLE 2" "not found - optional, the built-in Tcl engine is used instead"
fi

# Accelerators: NOT in the plugin zip. They arrive either as the per-OS
# binaries release asset (unpack it here as ./binaries/), or from a source
# build of the repository (./native/). Detect both.
BINDIR=""
for d in "$SRC/binaries" "$SRC"/vmdhole-binaries-*/ "$SRC/native"; do
    [ -x "$d/mole_tunnel_engine" ] && { BINDIR="${d%/}"; break; }
done
if [ -n "$BINDIR" ]; then
    note "tunnel engine" "found ($BINDIR)"
else
    note "tunnel engine" "not present - Tunnel mode falls back to Tcl (slower). Download the vmdhole-binaries asset for your OS and unpack it next to this script."
fi

# Trajectory data is NOT distributed with the plugin - see README.
note "trajectory data" "downloaded separately, see README (not bundled)"

if [ "$CHECK_ONLY" -eq 1 ]; then
    [ "$miss" -eq 0 ] && echo "All required components present." || echo "Required components are missing."
    exit "$miss"
fi
[ "$miss" -eq 0 ] || { echo; echo "Install aborted: VMD is required."; exit 1; }

echo
echo "Installing to $DEST/vmdhole"
mkdir -p "$DEST/vmdhole"
# Only what VMD loads: the plugin package, its licence/notice, and the binaries.
cp "$SRC/vmdhole/vmdhole.tcl"  "$DEST/vmdhole/"
cp "$SRC/vmdhole/pkgIndex.tcl" "$DEST/vmdhole/"
for f in NOTICE.md LICENSE-Apache-2.0.txt; do
    [ -f "$SRC/vmdhole/$f" ] && cp "$SRC/vmdhole/$f" "$DEST/vmdhole/"
done
echo "Installed."
if [ -n "$BINDIR" ]; then
    echo
    echo "Accelerator binaries found in $BINDIR:"
    echo "  point VMDHole at them under File > Settings (engine/sos_triangle paths),"
    echo "  or copy them next to your HOLE binaries so they are found automatically."
fi
echo
echo "Load it with:  vmd -e /dev/null   then  Extensions > Analysis > VMDHole"
echo "or add to ~/.vmdrc:  vmd_install_extension vmdhole vmdhole_tk \"Analysis/VMDHole\""
