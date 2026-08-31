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
# Forward the arguments: the underlying test takes a hole2/src checkout as $1
# and its SKIP message says to pass one - swallowing "$@" here made that advice
# a no-op and the group unrunnable through this wrapper.
exec "$T" "$@"
