# The paper's numbers

This directory is the quantitative record behind the VMDHole paper, and only
that. `benchmarks/results/` holds every reported CSV; each file's header
records exactly how, where and from which commit it was produced
(`env_manifest.txt` describes the machine).

The repository is the software. Everything that *generated* these numbers -
the benchmark harness (`reproduce.sh` and its stages), the comparison
fixtures, the case-study inputs and outputs, raw logs and rendered figures -
is archived in the companion data deposit so the record here stays small and
the replication kit stays complete. [DATA.md](DATA.md) lists the deposit's
contents and DOI, and where each trajectory used by the paper is published.

To replicate: download the deposit, unpack it over a checkout of this
repository, and run `paper/benchmarks/reproduce.sh` - the harness checks its
own preconditions and refuses to time anything on a broken or busy build.
