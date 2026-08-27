# MOLE 2: the tunnel similarity filter deletes the tunnel it selected

Reported against `github.com/sb-ncbr/MOLE` at `e0df21c` (25 Aug 2025), the tip
at the time of writing. Reproduced with Mono 6.12.0.199 on Linux.

This is why VMDHole's tunnel engine reports **8.44 Å** on 1MXT where MOLE 2
reports **12.04 Å** — the one documented deviation from MOLE's own output. The
deviation is pinned by `vmdhole/tests/mole_upstream_deviation.py` so it cannot
be "fixed" by accident.

## Summary

`TunnelComputation.FilterTunnels` decides which tunnel of a similar group to
keep, then removes the others with `List<Tunnel>.Remove`. That call uses
`EqualityComparer<Tunnel>.Default`, which resolves to `Tunnel`'s **value-based**
`Equals` (same lining residues). Where the filter needs object identity it gets
value equality, so the removal can delete the very tunnel the filter selected
and leave a look-alike in its place. The exported result is then a different
tunnel from the one the algorithm chose.

## Root cause

Two distinct problems. **The wrong survivor comes from the first; fixing the
second alone would not correct it.**

**1. Identity is required, value equality is used.** `FilterTunnels` reasons
about specific tunnel objects — it selects a pivot and marks particular others
for removal — but removes them with:

`WebChemistry.Tunnels.Core/TunnelComputation.cs:190`

```csharp
foreach (var t in toRemove) ret.Remove(t);   // ret is List<Tunnel>
```

`List<T>.Remove` performs a linear scan with `EqualityComparer<T>.Default`,
which for `Tunnel : IEquatable<Tunnel>` calls the value-based `Equals`. It
deletes the **first tunnel with a matching lining**, which need not be the
tunnel that was passed in.

**2. `Equals` and `GetHashCode` disagree.** `WebChemistry.Tunnels.Core/Tunnel.cs`

```csharp
public bool Equals(Tunnel other)            // VALUE semantics
{
    if (other == null) return false;
    if (this.lining.Count() != other.lining.Count()) return false;
    return Enumerable.Zip(this.lining, other.lining,
                          (l, r) => l.Identifier == r.Identifier).All(e => e);
}

public override int GetHashCode()           // REFERENCE semantics
{
    return base.GetHashCode();
}
```

This violates the documented contract that equal objects return equal hash
codes. Its practical effect here is on `HashSet<Tunnel> toRemove`: value-equal
tunnels are stored as separate entries, and the `toRemove.Contains(...)` guards
recognise only the same object.

## Observed

PDB 1MXT, default parameters, `<Origins Auto="1"/>`. One origin group of six
tunnels, five of which share a lining:

```
[0] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[1]  8.44/p13  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[2] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[3] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[4]  8.89/p13  lining(7) = 247 A,41 A,249 A,268 A,67 A,68 A,272 A
[5] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
```

The filter orders by path length and calls `RemoveLonger(shorter, ...)`, taking
8.44 as pivot and marking 8.89 and all four 12.04s. Instrumenting the removal
loop:

```
req  8.89 -> removes [4] 8.89   sameObject=True
req 12.04 -> removes [0] 12.04  sameObject=True
req 12.04 -> removes [0] 8.44   sameObject=False   <-- the pivot
req 12.04 -> removes [0] 12.04  sameObject=False
req 12.04 -> removes [0] 12.04  sameObject=False
survivors: 12.04
```

The policy is to keep the shorter tunnel of a similar pair; the export contains
the longer one. The same shape recurs on 1KX5 (cavity 18: two candidates,
identical 9-residue linings, MOLE prints the longer) — two structures, two
independent confirmations.

The outcome is stable across runs: `base.GetHashCode()` varies per run, the
output does not, because `HashSet<T>` enumerates its sequentially stored
entries in insertion order when nothing is removed. That matches the .NET
Framework reference source and the current .NET source, but is **not a public
guarantee** — read it as stable on the tested runtime, not as guaranteed.

## Minimal reproduction

`TunnelRemovalRepro.cs` in this directory, ~30 lines, no MOLE dependency:

    mcs -out:repro.exe TunnelRemovalRepro.cs && mono repro.exe

    filter decided to KEEP 8.44 and remove the other 5
    survivors: 12.04
    List.Contains(pivot) says      : True
    but is the pivot OBJECT there? : False
    BUG - the tunnel the filter chose to KEEP was deleted, a look-alike kept

Note the last two lines. `List.Contains(pivot)` reports true after the pivot has
been deleted, because the survivor shares its lining. **A test written with the
same value equality cannot detect this**; the assertion has to be
identity-based.

## Suggested fix

**Remove by identity.** `FilterTunnels` already assigns each tunnel a unique
per-call id at line 141 (`ret.ForEach((t, i) => t.Id = i);`), so:

```csharp
var toRemove = new HashSet<int>();
...
toRemove.Add(other.Id);
...
ret.RemoveAll(t => toRemove.Contains(t.Id));
```

This removes exactly the tunnels the filter selected, regardless of equality
semantics. **Separately, make `GetHashCode` honour `Equals`** — worth doing on
its own, but it does **not** fix the wrong survivor: `List.Remove` would still
match by value.

## Verified by rebuilding MOLE with the fix

The one-line identity-based removal was built and run, not proposed on paper:

    1MXT   stock  3.47 11.93 12.04 25.31
           fixed  3.47  8.44 11.93 25.31
    1KX5   stock  ... 7.72 7.84 8.75 ...   (60 tunnels)
           fixed  ... 7.61 7.72 8.75 ...   (60 tunnels)

On 1KX5 exactly one of the sixty tunnels moves; the other 59 are identical, so
the patch changes only the case it is meant to. **The fixed MOLE agrees with
VMDHole's engine exactly on both structures** — the strongest confirmation
available: if the cause were anything else, a one-line change to the removal
could not have closed the gap.

## Decisions taken in this repository

- The defect is **not replicated** in VMDHole's engine. Being wrong in the same
  way as the reference is a poor trade; VMDHole reports the tunnel MOLE's own
  filter selects.
- The 1MXT deviation is a **named regression** (`mole_upstream_deviation.py`)
  so it cannot silently disappear or grow.
- The acceptance bar for the port is "exact agreement with MOLE except
  documented upstream defects" — this is the one such defect.
