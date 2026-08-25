#!/bin/sh
# Thin wrapper so run_tests.sh picks up the HCAPEN cutoff-cache regression test,
# which lives next to the Fortran it tests (native/) rather than here.
# Both HCAPEN defects live in SAVE'd state ACROSS calls, so no whole-pipeline
# fixture can see them - verify.sh missed both.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
T="$DIR/../../native/connolly_patches/test_hcapen_cache.sh"
if [ ! -x "$T" ]; then
    echo "SKIP: native/connolly_patches/test_hcapen_cache.sh not found"
    exit 0
fi
exec "$T"
