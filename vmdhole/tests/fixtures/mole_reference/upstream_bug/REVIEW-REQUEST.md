# Review request: a suspected upstream bug in MOLE 2

I need a second pair of eyes before we act on this. Below is what we're building,
what I believe the bug is, and — most importantly — **what to check to prove me
wrong**. Three claims I made on the way to this conclusion turned out to be
false and were retracted (table in section 4), so please treat this as unproven
until you've run at least the first two checks.

---

## 1. What we're doing

We're porting **MOLE 2** (a tool that finds tunnels through proteins;
`github.com/sb-ncbr/MOLE`, C#/.NET) into a VMD plugin. Two implementations:
C for speed, Tcl so it works without the compiled binary. They produce
byte-identical output.

**The acceptance bar is exact agreement with MOLE's own output, except for
documented upstream defects.** Everything is validated against MOLE's own
exported files — profiles, lining layers, per-residue chemistry, cavities — and
we match to the last digit MOLE prints on every test structure. The suite is 24
checks, one of which pins the single documented deviation below.

**One structure disagrees**: myoglobin (PDB 1MXT). MOLE reports a tunnel of
length 12.04 Å; we report one of 8.44 Å. The other three tunnels agree exactly.
This document is about that one difference.

---

## 2. What I believe the bug is

### The code

`WebChemistry.Tunnels.Core/Tunnel.cs` — equality is **value-based**, hash is
**reference-based**:

```csharp
public bool Equals(Tunnel other)            // two tunnels are "equal" if they
{                                           // touch the same lining residues
    if (other == null) return false;
    if (this.lining.Count() != other.lining.Count()) return false;
    return Enumerable.Zip(this.lining, other.lining,
                          (l, r) => l.Identifier == r.Identifier).All(e => e);
}

public override int GetHashCode()
{
    return base.GetHashCode();              // identity hash — inconsistent
}
```

`WebChemistry.Tunnels.Core/TunnelComputation.cs:190`:

```csharp
foreach (var t in toRemove) ret.Remove(t);  // ret is List<Tunnel>
```

### Why that misbehaves

`List<T>.Remove` scans linearly using `EqualityComparer<T>.Default`, which for
`Tunnel : IEquatable<Tunnel>` calls the **value-based** `Equals`. So it deletes
the first tunnel with a *matching lining*, not the tunnel passed in.

Meanwhile `HashSet<Tunnel> toRemove` hashes first, so with a reference-based
hash it stores value-equal tunnels as separate entries, and in the observed run -
where the identity hashes all differ - the `toRemove.Contains(...)` guards a few
lines above recognise only the same *object*. A hash collision is a real, if
unlikely, exception to that.

**The wrong survivor comes specifically from using default VALUE equality where
`FilterTunnels` requires object identity.** The `Equals`/`GetHashCode` contract
violation is real and worth fixing on its own, but fixing `GetHashCode` alone
would NOT correct the wrong survivor - `List.Remove` would still match by value.

### The effect

The similarity filter selects a survivor; the removal loop deletes that
survivor and leaves a look-alike. On 1MXT, cavity 1, origin 3 — six tunnels,
five sharing a lining:

```
[0] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[1]  8.44/p13  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A     <- filter's pick
[2] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[3] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A
[4]  8.89/p13  lining(7) = 247 A,41 A,249 A,268 A,67 A,68 A,272 A
[5] 12.04/p18  lining(6) = 247 A,41 A,249 A,268 A,67 A,68 A

req  8.89 -> removes [4] 8.89   sameObject=True
req 12.04 -> removes [0] 12.04  sameObject=True
req 12.04 -> removes [0] 8.44   sameObject=False   <- deletes the filter's pick
req 12.04 -> removes [0] 12.04  sameObject=False
req 12.04 -> removes [0] 12.04  sameObject=False
survivors: 12.04
```

**Claim: our 8.44 Å is what MOLE's own filter selects. Their 12.04 Å is what
their removal loop leaves behind.**

---

## 3. What to check

Ordered by how much they'd change the conclusion.

### Check 1 — the standalone reproduction (2 minutes, no MOLE needed)

```sh
# mcs/mono are not on PATH; they live in the oracle env:
export PATH=<work>/mole_oracle_env/bin:$PATH
cd vmdhole/tests/fixtures/mole_reference/upstream_bug
mcs -out:/tmp/repro.exe TunnelRemovalRepro.cs && mono /tmp/repro.exe
```

Expected:

```
filter decided to KEEP 8.44 and remove the other 5
survivors: 12.04
List.Contains(pivot) says      : True
but is the pivot OBJECT there? : False
BUG - the tunnel the filter chose to KEEP was deleted, a look-alike kept
```

**What I want you to judge:** is my 30-line model actually faithful to MOLE's
code, or have I built a strawman that reproduces a bug MOLE doesn't have? The
model asserts three things — value-based `Equals`, reference-based
`GetHashCode`, and removal via `List.Remove`. Please confirm all three against
the real source rather than taking my quotes.

> Note the third and fourth output lines. My first version of this reproduction
> checked `List.Contains(pivot)` and printed **"OK, no bug"** — the check was
> fooled by the same broken equality that causes the bug. Worth remembering if
> you write your own test.

### Check 2 — is it upstream, or our local copy?

There are three copies of MOLE's source on this machine, and I initially read
the wrong one.

```sh
cd <work>/mole2_reference
git remote -v          # expect github.com/sb-ncbr/MOLE.git
git status --short     # expect empty (pristine)
git fetch && git log --oneline HEAD..origin/HEAD   # expect empty (at tip)
```

Then read `Tunnel.cs` and `TunnelComputation.cs:190` **in that clone**, not in
`mole_build/src` (which I modified with debug hooks) and not in
`webchemistry_reference`.

**What I want you to judge:** whether the code I'm quoting really is upstream's.

### Check 3 — how stable is it?

I originally said reproducing it would be risky because it depends on .NET
`HashSet` iteration order. That was overstated; the corrected position is that
the behaviour is stable on the tested runtime and expected for the
implementations inspected, but not guaranteed.

```sh
cd <work>/mole_build
for i in 1 2 3; do
  MOLE_REMOVE_DEBUG=1 ../mole_oracle_env/bin/mono mole2.exe tests/1MXT.xml 2>&1 \
    | grep -m1 "8.44/p13  hash"
  cut -d, -f2 out_1MXT/csv/tunnels.csv | tail -n +2 | tr '\n' ' '; echo
done
```

Expected: the hash differs every run, the output never does.

My reasoning: `HashSet<T>` enumerates its sequentially stored entries; `Add`
appends; nothing is ever removed from `toRemove`, so no freed slots are reused;
therefore enumeration follows insertion order. That matches the .NET Framework
reference source and the current .NET source.

Two caveats I do NOT want glossed over: enumeration order is not a public
guarantee, and an identity-hash collision would change `HashSet` membership and
hence the guards.

**What I want you to judge:** whether "stable on the tested runtime, expected
for the inspected implementations" is the right strength of claim, or whether
even that is too strong.

### Check 4 — could it be intentional?

The strongest counter-argument is that deleting every same-lining tunnel is
deliberate de-duplication. My reasons for rejecting that:

- the similarity filter has already computed which member to keep, and the
  removal deletes that member;
- 8.44 Å and 12.04 Å share lining residues but are 13- and 18-tetrahedron paths
  of very different length — same residues, different routes;
- a deliberate value-based dedup would need the `HashSet` guards to be
  value-based too, which the broken hash prevents;
- `GetHashCode` returning `base.GetHashCode()` under an empty `/// <summary>`
  looks like an unfinished override.

**What I want you to judge:** whether any of that is wishful reasoning toward
the answer that makes our port look right.

### Check 5 — reproduce it end to end (optional)

```sh
cd <work>/mole_build
../mole_oracle_env/bin/mono mole2.exe tests/1MXT.xml
cut -d, -f1,2 out_1MXT/csv/tunnels.csv      # MOLE: 11.93 12.04 25.31 3.47
```

Ours, same structure, same parameters:

```sh
cd <repo>
# atom table built from MOLE's own PDB, so both sides read the same input:
python3 vmdhole/tests/fixtures/mole_reference/extract_atoms_pdb.py \
        <work>/mole_build/tests/1MXT.pdb /tmp/1mxt.txt
./native/mole_tunnel_engine /tmp/1mxt.txt /tmp/out.txt
grep '^T ' /tmp/out.txt | awk '{printf "%.2f ", $4}'   # ours: 8.44 11.93 25.31 3.47
```

The rebuild recipe for MOLE (needed for the instrumented runs) is in
`vmdhole/NOTES/mole2-claims-audit.md`. Note `tc2.rsp`, not `tc.rsp`.

---

## 4. Where I've been wrong already

So you know where to press hardest. The verdict "this is a MOLE bug" has not
changed since I found it, but the reasoning around it has, three times:

| I said | Reality |
|---|---|
| MOLE's k-d tree doesn't return neighbours in distance order | False. `PriorityArray` is a sorted insert. Unrelated to this bug. |
| 1MXT is a tunnel *ordering* problem | False. The ordering was already correct; the tunnel *set* differs. |
| `Remove(8.89)` is the call that deletes 8.44 | False. 8.89 has a 7-residue lining, isn't equal to the others, and correctly removes only itself. It's the third `Remove(12.04)`. |
| Reproducing the bug would be machine-dependent | False, I believe — see Check 3. |

The first three were fixed by measuring instead of reading. The pattern I'd
watch for is me inferring a mechanism from source and stating it as a finding.

`vmdhole/NOTES/mole2-claims-audit.md` separates every "MOLE behaves like X"
claim in this port into *measured* (a mutation makes a comparison against MOLE
go red, or MOLE's own instrumented run prints it) versus *only read*. Three
behaviours in that file are still only read — no test structure triggers them.

---

## 5. Reviewed and decided

Independent review completed; the conclusion was accepted. Decisions taken:

- **Report upstream.** `REPORT.md` in this directory is the version to send. It
  now recommends removing by the unique per-call `Id` that `FilterTunnels`
  already assigns (`ret.ForEach((t, i) => t.Id = i);`, TunnelComputation.cs:141)
  - `HashSet<int>` plus `RemoveAll` by id - with the `GetHashCode` fix as a
  separate, and insufficient on its own, change.
- **Do NOT reproduce the defect** in our default implementation. Being wrong in
  the same way as the reference is a poor trade. If a genuine downstream
  zero-diff requirement appears, it gets an explicit compatibility mode rather
  than contaminating the correct default algorithm.
- **The 1MXT deviation is now a named regression** so it cannot be "fixed" by
  accident: `mole_upstream_deviation.py`, wired into the suite.
- **The acceptance bar is restated** as "exact agreement with MOLE except
  documented upstream defects", not "exact agreement".

## 6. The original open question (now answered above)

If you agree it's a bug:

1. **Do we replicate it** so we match MOLE exactly, or stay correct and document
   the difference? Replication is deterministic (Check 3), so it's tractable: we
   need each candidate's path lining at filter time, and a value-based delete
   instead of our current mark-and-skip. Roughly a day with tests.
2. **Do we report it upstream?** `REPORT.md` in this directory is written and
   ready to send, with the reproduction attached.

My inclination is to report it and *not* replicate — being wrong in the same way
as the reference is a poor trade for one tunnel on one structure — but that
depends on whether anyone downstream diffs our output against MOLE's and expects
zero differences.
