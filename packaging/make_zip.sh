#!/bin/sh
# Builds packaging/vmdhole.zip straight from the source trees - no hand-maintained
# duplicate directory in between.
#
# There used to be a whole second copy of the plugin, packaging/vmdhole/, kept in
# sync by hand on every change and separately zipped. It existed only because
# this script once zipped whatever sat in that directory rather than building
# from source, and CI's own release job just re-uploaded that checked-in zip
# without ever rebuilding it - so "keep the duplicate in sync" was a real,
# recurring cost with no reproducibility benefit (a fresh build from source
# proves the package matches HEAD; a hand-copied one only proves someone
# remembered to copy). This script now assembles the same content into a
# throwaway staging directory and zips that, so there is exactly one source of
# truth: vmdhole/, native/, docs/, and the top-level
# install.sh/README/LICENSE.
#
# The archive is what a USER installs, and ONLY that: the plugin package
# (vmdhole/ - the script, pkgIndex, licences, and the two small tutorial
# structures), install.sh, README and LICENSE. Nothing else:
#   - compiled accelerators are the per-OS release assets built by
#     .github/workflows/binaries.yml (unpack one next to install.sh as
#     binaries/ and the installer wires it up; the plugin runs without them
#     on its built-in Tcl engines, only slower, and Tunnel mode needs the
#     engine);
#   - the C/Fortran SOURCE, docs source, tests and benchmarks live in the
#     repository - a user who wants to build or develop clones it;
#   - the documentation is hosted (see README) and not duplicated here.
# Re-packing the repository into every user download was the old behaviour;
# the package is now the product, not the project.
#
# Built from `git ls-files`, not a filesystem `find`: the working tree also
# holds untracked local build directories (native/stock_build/,
# sphproc_build/, build/ - compiled HOLE checkouts and binaries, tens of MB)
# that a raw `find` swept straight into the archive - one early version of
# this script produced a 353 MB, 5351-file zip. git ls-files only ever returns
# what is actually in the repository.
set -e
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO=$(CDPATH= cd -- "$DIR/.." && pwd)
cd "$REPO"

STAGE="$DIR/build/VMDHole"
rm -rf "$DIR/build"
mkdir -p "$STAGE"

copy_tracked() {
    # copy_tracked <git-ls-files-prefix> [exclude-prefix ...]
    #
    # Skips a tracked path that is not actually on disk, rather than failing:
    # `git ls-files` lists what HEAD (plus the index) has, which can include a
    # file an uncommitted local change has already deleted. That is a real,
    # separate, in-progress state this script has no business resolving - it
    # only decides what belongs in the ARCHIVE, not what belongs in a commit.
    prefix=$1; shift
    git ls-files "$prefix" | while IFS= read -r f; do
        skip=0
        for ex in "$@"; do
            # A trailing "/" excludes a whole directory (prefix match); anything
            # else is an EXACT file match, or "native/hydro_project" would
            # also swallow "native/hydration/hydro_project.c" as a false prefix hit -
            # measured: it did, on the first version of this exclusion list.
            case "$ex" in
                */) case "$f" in "$ex"*) skip=1; break;; esac ;;
                *)  [ "$f" = "$ex" ] && { skip=1; break; } ;;
            esac
        done
        [ "$skip" = 1 ] && continue
        [ -f "$f" ] || { echo "  (skipping $f - tracked but not on disk)" >&2; continue; }
        mkdir -p "$STAGE/$(dirname "$f")"
        cp "$f" "$STAGE/$f"
    done
}

copy_tracked vmdhole/ vmdhole/NOTES/ vmdhole/tests/ vmdhole/hole_tcl/ vmdhole/HANDOFF.md

# Top level: licence, readme, installer - not under vmdhole/, so
# copy_tracked's prefix match does not reach them. (The logo is not shipped:
# it exists for the repository front page.)
cp LICENSE "$STAGE/LICENSE"
cp README.md "$STAGE/README.md"
cp install.sh "$STAGE/install.sh"
chmod +x "$STAGE/install.sh"

# Optional: VMDHOLE_BINARIES_DIR names a directory of compiled accelerators
# (the per-OS "wheel" set binaries.yml builds at -O0). When given, it is
# bundled as vmdhole/binaries/ - the exact layout install.sh detects - so CI
# can assemble ONE download per OS: plugin + matching binaries. Unset (the
# default, and what the release-integrity test checks) the zip stays the
# portable plugin-only package.
if [ -n "${VMDHOLE_BINARIES_DIR:-}" ]; then
    [ -d "$VMDHOLE_BINARIES_DIR" ] || { echo "VMDHOLE_BINARIES_DIR is not a directory: $VMDHOLE_BINARIES_DIR" >&2; exit 1; }
    mkdir -p "$STAGE/binaries"
    cp "$VMDHOLE_BINARIES_DIR"/* "$STAGE/binaries/"
    chmod +x "$STAGE/binaries"/* 2>/dev/null || true
fi

rm -f "$DIR/vmdhole.zip"
( cd "$DIR/build" && find VMDHole -type f | LC_ALL=C sort | zip -q -X "$DIR/vmdhole.zip" -@ )
echo "vmdhole.zip: $(unzip -l "$DIR/vmdhole.zip" | tail -1 | awk '{print $2}') files, $(du -h "$DIR/vmdhole.zip" | cut -f1)"
