// constrained_mcmc.cpp -- degree-conditional nulls of Definition S1 (a), (c), (d)
//
// Markov chain on directed binary networks (no self-loops, no multiple links)
// with fixed in- and out-degrees.  Proposals (all symmetric):
//   0  double-edge swap      (a->b, c->d)          -> (a->d, c->b)
//   1  reciprocal-pair swap  (a<->b, c<->d)        -> (a<->d, c<->b)
//   2  reversal of a directed three-cycle with no reciprocal link
// mode 0: no further constraint (exact rewiring, Definition S1(a); proposals 0 and 2)
// mode 1: number of reciprocal pairs fixed (Definition S1(c))
// mode 2: number of reciprocal partners of every node fixed (Definition S1(d))
// A proposal that would violate the constraint is rejected, so the stationary
// distribution is uniform on the part of the constrained set reachable from the
// starting network.  Random numbers come from R's generator (set.seed in R).
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

namespace {

struct Net {
  int N, E;
  std::vector<unsigned char> A;   // A[i * N + j] = 1 for a link i -> j
  std::vector<int> src, dst, eid; // edge list and index of each link
  inline unsigned char a(int i, int j) const { return A[i * N + j]; }
  inline bool rec(int i, int j) const { return A[i * N + j] && A[j * N + i]; }
  void move_edge(int e, int i, int j) {          // edge e becomes i -> j
    int oi = src[e], oj = dst[e];
    A[oi * N + oj] = 0; eid[oi * N + oj] = -1;
    src[e] = i; dst[e] = j; A[i * N + j] = 1; eid[i * N + j] = e;
  }
};

inline int rint(int n) { int k = (int) (unif_rand() * n); return k < n ? k : n - 1; }

// double-edge swap with the constraint check of the chosen mode
bool prop_swap(Net& g, int mode) {
  int e1 = rint(g.E), e2 = rint(g.E);
  if (e1 == e2) return false;
  int a = g.src[e1], b = g.dst[e1], c = g.src[e2], d = g.dst[e2];
  if (a == c || b == d || a == d || c == b) return false;
  if (g.a(a, d) || g.a(c, b)) return false;
  int px[4] = {a, c, a, c}, py[4] = {b, d, d, b};
  int before[4], after[4];
  if (mode > 0) for (int k = 0; k < 4; k++) before[k] = g.rec(px[k], py[k]);
  g.move_edge(e1, a, d); g.move_edge(e2, c, b);
  if (mode == 0) return true;
  for (int k = 0; k < 4; k++) after[k] = g.rec(px[k], py[k]);
  bool ok = true;
  if (mode == 1) {
    int s = 0; for (int k = 0; k < 4; k++) s += after[k] - before[k];
    ok = (s == 0);
  } else {
    int nodes[4] = {a, b, c, d};
    for (int v = 0; v < 4 && ok; v++) {
      int s = 0;
      for (int k = 0; k < 4; k++)
        if (px[k] == nodes[v] || py[k] == nodes[v]) s += after[k] - before[k];
      ok = (s == 0);
    }
  }
  if (!ok) { g.move_edge(e1, a, b); g.move_edge(e2, c, d); }
  return ok;
}

// swap of two reciprocal pairs: preserves degrees, reciprocal pairs and partners
bool prop_recswap(Net& g) {
  int e1 = rint(g.E), e2 = rint(g.E);
  if (e1 == e2) return false;
  int a = g.src[e1], b = g.dst[e1], c = g.src[e2], d = g.dst[e2];
  if (a == c || a == d || b == c || b == d) return false;
  if (!g.a(b, a) || !g.a(d, c)) return false;
  if (g.a(a, d) || g.a(d, a) || g.a(c, b) || g.a(b, c)) return false;
  int f1 = g.eid[b * g.N + a], f2 = g.eid[d * g.N + c];
  g.move_edge(e1, a, d); g.move_edge(f1, d, a);
  g.move_edge(e2, c, b); g.move_edge(f2, b, c);
  return true;
}

// reversal of a three-cycle a->b->c->a without reciprocal links
bool prop_cycle(Net& g) {
  int e1 = rint(g.E), c = rint(g.N);
  int a = g.src[e1], b = g.dst[e1];
  if (c == a || c == b) return false;
  if (!g.a(b, c) || !g.a(c, a)) return false;
  if (g.a(b, a) || g.a(c, b) || g.a(a, c)) return false;
  int e2 = g.eid[b * g.N + c], e3 = g.eid[c * g.N + a];
  g.move_edge(e1, b, a); g.move_edge(e2, c, b); g.move_edge(e3, a, c);
  return true;
}

bool step(Net& g, int mode) {
  double u = unif_rand();
  if (mode == 0) return (u < 0.9) ? prop_swap(g, 0) : prop_cycle(g);
  if (u < 0.45) return prop_swap(g, mode);
  if (u < 0.90) return prop_recswap(g);
  return prop_cycle(g);
}

} // namespace

// A0[i, j] = 1 for a link i -> j (the direction convention does not matter for
// the statistics used here, which depend on A only through A * t(A)).
// [[Rcpp::export]]
List constrained_chain_cpp(IntegerMatrix A0, NumericVector w, int mode, int n_draws,
                           double burn_mult = 500, double thin_mult = 100,
                           bool keep_mean = false, bool keep_codes = false) {
  int N = A0.nrow();
  Net g; g.N = N; g.A.assign(N * N, 0); g.eid.assign(N * N, -1);
  for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
    if (i != j && A0(i, j) != 0) {
      g.A[i * N + j] = 1; g.eid[i * N + j] = g.src.size();
      g.src.push_back(i); g.dst.push_back(j);
    }
  g.E = g.src.size();
  double hhi = 0; for (int i = 0; i < N; i++) hhi += w[i] * w[i];
  NumericVector psi(n_draws); IntegerVector npairs(n_draws), codes(keep_codes ? n_draws : 0);
  NumericMatrix mean_rec(keep_mean ? N : 0, keep_mean ? N : 0);
  long long burn = (long long) (burn_mult * g.E), thin = (long long) (thin_mult * g.E);
  if (thin < 1) thin = 1;
  long long acc = 0, tot = 0;
  if (g.E >= 2) for (long long s = 0; s < burn; s++) { acc += step(g, mode); tot++; }
  for (int k = 0; k < n_draws; k++) {
    if (g.E >= 2) for (long long s = 0; s < thin; s++) { acc += step(g, mode); tot++; }
    double sw = 0; int np = 0;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++)
      if (g.rec(i, j)) { sw += 2 * w[i] * w[j]; np++;
        if (keep_mean) { mean_rec(i, j) += 1; mean_rec(j, i) += 1; } }
    psi[k] = sw / (1 - hhi); npairs[k] = np;
    if (keep_codes) { int code = 0, bit = 0;
      for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) {
        if (g.a(i, j)) code |= (1 << bit); bit++; }
      codes[k] = code; }
  }
  if (keep_mean) mean_rec = mean_rec / (double) n_draws;
  // final degrees, as a check that the chain preserved them
  IntegerVector kout(N), kin(N);
  for (int e = 0; e < g.E; e++) { kout[g.src[e]]++; kin[g.dst[e]]++; }
  return List::create(_["psi"] = psi, _["n_pairs"] = npairs, _["mean_rec"] = mean_rec,
                      _["codes"] = codes, _["accept"] = tot > 0 ? (double) acc / tot : NA_REAL,
                      _["k_out"] = kout, _["k_in"] = kin, _["n_links"] = g.E);
}
