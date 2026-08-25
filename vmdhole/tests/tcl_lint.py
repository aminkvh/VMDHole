#!/usr/bin/env python3
"""Static checks for Tcl defects that SOURCE CLEANLY but throw when called.

Every check here exists because the class actually bit this project, or is one
brace away from it. `info complete` and a successful `source` prove nothing
about any of them - the file loads fine and the proc explodes the first time a
user reaches it.

  switch-comment   a switch body is parsed as a LIST, so a '#' line between its
                   cases is not a comment: its words become spurious pattern/
                   body pairs, and a quoted phrase followed by a comma is
                   invalid list syntax outright. This shipped in property_meta
                   and made it throw for EVERY property until a user hit it.
                   Comments inside a case's braced body are legal and common,
                   so only top-level lines between cases are flagged.
  switch-parity    a switch body with an odd element count throws "extra switch
                   pattern with no body" on the first call.
  dead-callback    -command / bind naming a ::VMDHole:: proc that does not
                   exist. The widget builds fine and does nothing when clicked.
  undeclared-var   a namespace variable read inside a proc without `variable`.
                   Reads as empty (or throws), silently. Hit for real when the
                   HOLE panel moved into a notebook page and two procs kept
                   using $_runpanel without declaring it.

Usage:  tcl_lint.py FILE      -> prints findings, exit 1 if any
"""
import re
import sys


def proc_bodies(lines):
    """(start_line_index, name, [body lines]) with brace-matched extent."""
    out = []
    for i, l in enumerate(lines):
        m = re.match(r'\s*proc\s+(\S+)\s*\{', l)
        if not m:
            continue
        depth = 0
        started = False
        body = []
        for j in range(i, len(lines)):
            body.append(lines[j])
            depth += lines[j].count("{") - lines[j].count("}")
            if "{" in lines[j]:
                started = True
            if started and depth <= 0:
                break
        out.append((i, m.group(1), body))
    return out


def switch_bodies(lines):
    """(switch_line_no, [(line_no, text)]) for each switch whose body opens on
    the same line. Yields only lines at the body's TOP level."""
    out = []
    i = 0
    while i < len(lines):
        l = lines[i]
        if re.search(r'(?<![\w-])switch\b', l) and l.rstrip().endswith("{"):
            base = l.count("{") - l.count("}")
            depth = base
            j = i + 1
            top, whole = [], []
            while j < len(lines) and depth > 0:
                if depth == base:
                    top.append((j + 1, lines[j]))
                whole.append(lines[j])
                depth += lines[j].count("{") - lines[j].count("}")
                j += 1
            out.append((i + 1, top, whole))
            i = j
            continue
        i += 1
    return out


def main():
    path = sys.argv[1]
    lines = open(path, encoding="utf-8").read().split("\n")
    src = "\n".join(lines)
    findings = []

    procs = {name: i + 1 for i, name, _ in proc_bodies(lines)}

    # --- switch-comment + switch-parity -------------------------------------
    for sw_ln, top, whole in switch_bodies(lines):
        for ln, txt in top:
            if re.match(r'^\s*#', txt):
                findings.append(("switch-comment", ln,
                                 f"comment between cases of the switch at L{sw_ln}: "
                                 f"{txt.strip()[:60]}"))
        # parity: only meaningful if the body is list-parseable at all, and the
        # comment check above already reports the un-parseable ones
        if not any(re.match(r'^\s*#', t) for _, t in top):
            body = "\n".join(whole)
            if body.count("{") == body.count("}"):
                toks = re.findall(r'\{|\}', body)
                depth = 0
                elems = 0
                for ln, txt in top:
                    s = txt.strip()
                    if not s:
                        continue
                    elems += len(re.findall(r'^\S+', s))
                # cheap parity proxy: count top-level case labels vs bodies
                labels = sum(1 for _, t in top if re.match(r'^\s*\S+\s*\{', t))
                bodies = sum(t.count("{") for _, t in top)
                if labels and bodies and labels != bodies:
                    findings.append(("switch-parity", sw_ln,
                                     f"{labels} case labels vs {bodies} bodies - "
                                     f"check for a missing/extra brace"))

    # --- dead-callback -------------------------------------------------------
    for i, l in enumerate(lines):
        for m in re.finditer(r'-command\s+(?:\{\s*)?(::VMDHole::[A-Za-z_][\w:]*)', l):
            if m.group(1) not in procs:
                findings.append(("dead-callback", i + 1,
                                 f"-command {m.group(1)} - no such proc"))
        for m in re.finditer(r'bind\s+\S+\s+<[^>]+>\s+(::VMDHole::[A-Za-z_][\w:]*)', l):
            if m.group(1) not in procs:
                findings.append(("dead-callback", i + 1,
                                 f"bind {m.group(1)} - no such proc"))

    # --- undeclared-var ------------------------------------------------------
    nsvars = set()
    for stmt in re.split(r'[;\n]', src):
        m = re.match(r'\s*variable\s+([A-Za-z_]\w*)', stmt)
        if m:
            nsvars.add(m.group(1))

    for idx, name, body in proc_bodies(lines):
        code = "\n".join(l for l in body if not l.strip().startswith("#"))
        decl = set()
        for stmt in re.split(r'[;\n]', code):       # `variable a; variable b`
            m = re.match(r'\s*variable\s+([A-Za-z_]\w*)', stmt)
            if m:
                decl.add(m.group(1))
        pm = re.match(r'\s*proc\s+\S+\s*\{([^}]*)\}', lines[idx])
        params = set(re.findall(r'[A-Za-z_]\w*', pm.group(1))) if pm else set()
        local = set()
        # `set x`, `[set x`, and `for {set x 0}` - the brace form matters
        local |= set(re.findall(r'[\s\[{]set\s+([A-Za-z_]\w*)', code))
        local |= set(re.findall(r'^\s*set\s+([A-Za-z_]\w*)', code, re.M))
        local |= set(w for g in re.findall(r'foreach\s+\{([^}]*)\}', code)
                     for w in re.findall(r'[A-Za-z_]\w*', g))
        local |= set(re.findall(r'foreach\s+([A-Za-z_]\w*)', code))
        local |= set(w for g in re.findall(r'lassign\s+\S+\s+([^\n\]]*)', code)
                     for w in re.findall(r'[A-Za-z_]\w*', g))
        local |= set(re.findall(r'upvar\s+\S+\s+([A-Za-z_]\w*)', code))
        local |= set(w for g in re.findall(r'dict\s+for\s+\{([^}]*)\}', code)
                     for w in re.findall(r'[A-Za-z_]\w*', g))
        for v in sorted(set(re.findall(r'\$([A-Za-z_]\w*)', code))):
            if v in nsvars and v not in decl and v not in params and v not in local:
                findings.append(("undeclared-var", idx + 1,
                                 f"{name} reads ${v} without `variable {v}`"))

    # no-BS audit: narrative/history comments. Ratcheted - the existing backlog
    # is allowed, growth is not.
    nobs = check_nobs(path, lines)
    if len(nobs) > NOBS_BASELINE:
        for ln, txt in nobs[-(len(nobs) - NOBS_BASELINE):]:
            findings.append(("no-BS", ln, f"narrative/history comment: {txt}"))

    if not findings:
        print("  clean: no switch-comment, switch-parity, dead-callback, "
              f"undeclared-var or no-BS findings ({len(nobs)}/{NOBS_BASELINE} "
              "narrative comments, not growing)")
        return 0
    by_kind = {}
    for kind, ln, msg in findings:
        by_kind.setdefault(kind, []).append((ln, msg))
    for kind in sorted(by_kind):
        print(f"  {kind}: {len(by_kind[kind])}")
        for ln, msg in by_kind[kind][:12]:
            print(f"      L{ln}: {msg}")
    return 1


# ---------------------------------------------------------------- no-BS audit
# Comments state what the code does and why, not what it used to do. Narrative,
# development history and re-telling a fixed bug are noise for every future
# reader. Baselined: the count may not GROW.
NOBS_PATTERNS = [
    # what the code used to do
    r"\bthis used to\b",
    r"\bit used to\b",
    r"\bused to (?:be|say|take|scan|relax|put|fall|return|gate|reposition|hardcode|force|claim|stop|rewrite|ask|clip|go|auto-approve|make|surface|come)\b",
    r"\bwas written when\b",
    r"\bthe exact failure\b",
    r"\bwhat this was\b",
    r"\blooked dead\b",
    r"\bmade every (?:global|inherited|trajectory)\b",
    r"\bCORRECTED, was wrong\b",
    r"\bSTALE PARAGRAPH\b",
    r"\bpreviously\b",
    r"\bback when\b",
    # session archaeology: which task, which day, which planning doc.
    # A reader of the shipped file has none of those.
    r"\b(?:task|defect batch|batch item)\s*#?\s*\d",
    r"\bitem\s+\d+\s+(?:already|fixed|added|is)\b",
    r"\b20\d\d-\d\d-\d\d\b",
    r"NOTES/",
    # talking to a colleague rather than describing the code
    r"(?<![A-Za-z])(?:we|I|my|our)(?![A-Za-z'])",
    r"\bnote that\b",
    r"\bkeep in mind\b",
    r"\bbeware\b",
]
# Tree count when each rule landed; may not GROW. The first ten are the original
# narrative rule; the rest were added by the shipped-code audit, which cleared
# every task/date/NOTES reference and left the first-person backlog measured.
NOBS_BASELINE = 166


def check_nobs(path, lines):
    import re as _re
    pat = _re.compile("|".join(NOBS_PATTERNS), _re.I)
    hits = []
    for i, l in enumerate(lines, 1):
        t = l.strip()
        if not t.startswith("#"):
            continue
        if pat.search(t):
            hits.append((i, t[:110]))
    return hits


if __name__ == "__main__":
    sys.exit(main())
