using System; using System.Collections.Generic; using System.Linq;
// Minimal reproduction of MOLE 2's tunnel-removal defect.
// Mirrors Tunnel: Equals compares the lining residue list (VALUE semantics),
// GetHashCode is left as base.GetHashCode() (REFERENCE semantics).
class Tunnel : IEquatable<Tunnel> {
    public string Lining; public double Length;
    public Tunnel(string lining, double len) { Lining = lining; Length = len; }
    public bool Equals(Tunnel o) { return o != null && Lining == o.Lining; }
    public override bool Equals(object o) { return o is Tunnel && Equals((Tunnel)o); }
    public override int GetHashCode() { return base.GetHashCode(); }   // <-- the defect
    public override string ToString() { return Length.ToString("0.00"); }
}
class Repro {
  static void Main() {
    const string L1 = "247A,41A,249A,268A,67A,68A";       // shared by five tunnels
    const string L2 = "247A,41A,249A,268A,67A,68A,272A";  // 8.89 has one residue more
    var ret = new List<Tunnel> {                      // 1MXT, cavity 1, origin 3
        new Tunnel(L1,12.04), new Tunnel(L1,8.44), new Tunnel(L1,12.04),
        new Tunnel(L1,12.04), new Tunnel(L2,8.89),  new Tunnel(L1,12.04) };
    // The similarity filter picks 8.44 as the pivot and marks the rest.
    var pivot = ret[1];
    var toRemove = new HashSet<Tunnel>();
    foreach (var t in ret) if (!ReferenceEquals(t, pivot)) toRemove.Add(t);
    Console.WriteLine("filter decided to KEEP {0} and remove the other {1}",
                      pivot, toRemove.Count);
    foreach (var t in toRemove) ret.Remove(t);        // TunnelComputation.cs:190
    Console.WriteLine("survivors: {0}",
        string.Join(", ", ret.Select(t => t.ToString()).ToArray()));
    // NOTE: ret.Contains(pivot) reports TRUE here even though the pivot object
    // is gone - the surviving 12.04 has the same lining, so the broken equality
    // fools the check as well. Identity is the only reliable test.
    Console.WriteLine("List.Contains(pivot) says      : {0}", ret.Contains(pivot));
    Console.WriteLine("but is the pivot OBJECT there? : {0}",
                      ret.Any(t => ReferenceEquals(t, pivot)));
    Console.WriteLine(ret.Any(t => ReferenceEquals(t, pivot))
        ? "OK - the kept tunnel survived"
        : "BUG - the tunnel the filter chose to KEEP was deleted, a look-alike kept");
  }
}
