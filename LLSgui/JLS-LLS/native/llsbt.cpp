// llsbt: a lifesrc/JLS-style backtracking search for LLS problems.
//
// Cells are decided in a fixed order (dead first) with chronological
// backtracking. Each cell transition is propagated exactly through a
// precomputed table over all 3^10 partial assignments of
// (predecessor, 8 neighbours, cell), so any value forced by a partly known
// neighbourhood is found immediately; this works for any isotropic rule.
// All of LLS's other constraints (symmetry links, population, change limits,
// ...) arrive as ordinary clauses and are propagated with watched literals.
//
// Usage: llsbt <cnf> <watch-file> <transitions-file> [options]
//   cnf:              clauses other than the transition rule
//   watch-file:       line 1 = variables to stream (in order),
//                     line 2 = variables that distinguish solutions
//   transitions-file: line 1 = 512 rule bits (B for masks 0..255, then S),
//                     line 2 = decision order (literals to try first),
//                     then one transition per line: pred n0..n7 cell,
//                     each a DIMACS literal, or T / F for constants
//   --interval MS     snapshot interval (default 250)
//   --solutions N     number of solutions (0 = all, default 1)
//   --threads N       search with N threads (the tree is split into parts)
//   --probe           before searching, try each cell both ways and fix the values
//                     that fail (JLS does this while preparing a search)
//
// stdin commands and stdout protocol are the same as llslive's.

#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

using Clock = std::chrono::steady_clock;

static std::atomic<bool> g_stop{false};
static std::atomic<bool> g_pause{false};
static std::atomic<bool> g_snap{false};
static std::mutex g_out;
static Clock::time_point g_start;

static double elapsed() {
  return std::chrono::duration<double>(Clock::now() - g_start).count();
}

static void emit(const std::string &s) {
  std::lock_guard<std::mutex> lock(g_out);
  fputs(s.c_str(), stdout);
  fputc('\n', stdout);
  fflush(stdout);
}

static const int TRUE_LIT = 1 << 30; // constant true in a transition
static const int FALSE_LIT = 0;      // constant false in a transition
static const int NPOS = 10;          // pred, 8 neighbours, cell
static const int TABLE_SIZE = 59049; // 3^10

struct Entry {
  uint8_t ok;           // some completion obeys the rule
  uint16_t force_true;  // positions true in every completion (unknown ones)
  uint16_t force_false; // positions false in every completion
};

struct Search {
  int nvars = 0;
  std::vector<std::vector<int>> clauses;
  std::vector<std::vector<int>> watches; // per literal index -> clause ids
  std::vector<std::array<int, NPOS>> trans;
  std::vector<std::vector<int>> occurs; // var -> transition ids
  std::vector<Entry> table;
  int pow3[NPOS];

  std::vector<int8_t> value; // per var: 1 true, -1 false, 0 unassigned
  std::vector<int> level;
  std::vector<int> trail;
  std::vector<int> trail_lim;     // trail size at each decision
  std::vector<char> flipped;      // per decision level
  size_t qhead = 0;

  std::vector<int> order; // decision literals, tried in this order
  size_t order_hint = 0;  // all order entries before this are assigned
  std::vector<char> is_block;

  long long decisions = 0, backtracks = 0;

  static int lit_index(int lit) { return 2 * abs(lit) + (lit < 0); }
  int lit_value(int lit) const { // 1 true, -1 false, 0 unknown
    int v = value[abs(lit)];
    return lit > 0 ? v : -v;
  }

  void build_table(const std::vector<int> &rule_bits) {
    pow3[0] = 1;
    for (int i = 1; i < NPOS; i++) pow3[i] = pow3[i - 1] * 3;
    // valid[full] for the 1024 total assignments: bit 0 = pred, 1..8 = neighbours, 9 = cell
    std::vector<char> valid(1024);
    for (int a = 0; a < 1024; a++) {
      int pred = a & 1, nbrs = (a >> 1) & 255, cell = (a >> 9) & 1;
      valid[a] = rule_bits[pred * 256 + nbrs] == cell;
    }
    table.assign(TABLE_SIZE, Entry{0, 0, 0});
    for (int idx = 0; idx < TABLE_SIZE; idx++) {
      int known_mask = 0, known_bits = 0, unknown[NPOS], nu = 0;
      for (int i = 0, r = idx; i < NPOS; i++, r /= 3) {
        int s = r % 3; // 0 unknown, 1 false, 2 true
        if (s == 0) unknown[nu++] = i;
        else {
          known_mask |= 1 << i;
          if (s == 2) known_bits |= 1 << i;
        }
      }
      uint16_t seen_true = 0, seen_false = 0;
      bool ok = false;
      for (int c = 0; c < (1 << nu); c++) {
        int a = known_bits;
        for (int j = 0; j < nu; j++)
          if (c >> j & 1) a |= 1 << unknown[j];
        if (!valid[a]) continue;
        ok = true;
        for (int j = 0; j < nu; j++) {
          if (c >> j & 1) seen_true |= 1 << unknown[j];
          else seen_false |= 1 << unknown[j];
        }
      }
      Entry &e = table[idx];
      e.ok = ok;
      if (ok) {
        int unknown_mask = ((1 << NPOS) - 1) & ~known_mask;
        e.force_true = unknown_mask & seen_true & ~seen_false;
        e.force_false = unknown_mask & seen_false & ~seen_true;
      }
    }
  }

  // Assign a literal true. Returns false on conflict.
  bool assign(int lit) {
    int v = abs(lit);
    int cur = lit_value(lit);
    if (cur == 1) return true;
    if (cur == -1) return false;
    value[v] = lit > 0 ? 1 : -1;
    level[v] = (int) trail_lim.size();
    trail.push_back(lit);
    return true;
  }

  bool check_transition(int id) {
    const auto &t = trans[id];
    int idx = 0;
    for (int i = 0; i < NPOS; i++) {
      int l = t[i], s;
      if (l == FALSE_LIT) s = 1;
      else if (l == TRUE_LIT) s = 2;
      else {
        int v = lit_value(l);
        s = v == 0 ? 0 : v > 0 ? 2 : 1;
      }
      idx += s * pow3[i];
    }
    const Entry &e = table[idx];
    if (!e.ok) return false;
    if (e.force_true | e.force_false) {
      for (int i = 0; i < NPOS; i++) {
        if (e.force_true >> i & 1) { if (!assign(t[i])) return false; }
        else if (e.force_false >> i & 1) { if (!assign(-t[i])) return false; }
      }
    }
    return true;
  }

  // Unit propagation over clauses and transitions. Returns false on conflict.
  bool propagate() {
    while (qhead < trail.size()) {
      int lit = trail[qhead++];
      // clauses watching the now-false literal -lit
      auto &ws = watches[lit_index(-lit)];
      for (size_t k = 0; k < ws.size();) {
        int cid = ws[k];
        auto &c = clauses[cid];
        if (c[0] == -lit) std::swap(c[0], c[1]);
        // c[1] is the false watch
        if (lit_value(c[0]) == 1) { k++; continue; }
        bool moved = false;
        for (size_t j = 2; j < c.size(); j++) {
          if (lit_value(c[j]) != -1) {
            std::swap(c[1], c[j]);
            watches[lit_index(c[1])].push_back(cid);
            ws[k] = ws.back();
            ws.pop_back();
            moved = true;
            break;
          }
        }
        if (moved) continue;
        if (!assign(c[0])) return false; // all false: conflict
        k++;
      }
      for (int id : occurs[abs(lit)])
        if (!check_transition(id)) return false;
    }
    return true;
  }

  void undo_to(size_t trail_size) {
    while (trail.size() > trail_size) {
      int v = abs(trail.back());
      value[v] = 0;
      trail.pop_back();
      // order_hint must not skip variables that became unassigned
    }
    qhead = trail.size();
    order_hint = 0;
  }

  // Chronological backtracking: undo levels until one can be flipped.
  // Returns false when the search space is exhausted.
  bool backtrack() {
    backtracks++;
    while (!trail_lim.empty()) {
      size_t start = trail_lim.back();
      int decision = trail[start];
      bool was_flipped = flipped.back();
      undo_to(start);
      trail_lim.pop_back();
      flipped.pop_back();
      if (!was_flipped) {
        trail_lim.push_back(trail.size());
        flipped.push_back(1);
        assign(-decision);
        if (propagate()) return true;
        // conflict right after flipping: keep unwinding
        backtracks++;
      }
    }
    return false;
  }

  int next_decision() {
    while (order_hint < order.size() && value[abs(order[order_hint])] != 0) order_hint++;
    if (order_hint < order.size()) return order[order_hint];
    return 0;
  }

  // After a solution: drop decision levels on variables that don't distinguish
  // solutions, so the next solution differs in a distinguishing variable.
  void skip_equivalent_solutions() {
    while (!trail_lim.empty() && !is_block[abs(trail[trail_lim.back()])]) {
      undo_to(trail_lim.back());
      trail_lim.pop_back();
      flipped.pop_back();
    }
  }
};

static std::vector<int> parse_ints(const std::string &line) {
  std::vector<int> out;
  std::istringstream in(line);
  int x;
  while (in >> x) out.push_back(x);
  return out;
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: llsbt <cnf> <watch-file> <transitions-file> [options]\n");
    return 1;
  }
  g_start = Clock::now();
  int interval_ms = 250;
  long long max_solutions = 1;
  bool probe = false;
  int threads = 1;
  for (int i = 4; i < argc; i++) {
    std::string a = argv[i];
    if (a == "--interval" && i + 1 < argc) interval_ms = atoi(argv[++i]);
    else if (a == "--solutions" && i + 1 < argc) max_solutions = atoll(argv[++i]);
    else if (a == "--probe") probe = true;
    else if (a == "--threads" && i + 1 < argc) threads = std::max(1, atoi(argv[++i]));
    else if ((a == "--seed" || a == "--opt") && i + 1 < argc) i++; // accepted for compatibility
    else {
      fprintf(stderr, "unknown argument %s\n", argv[i]);
      return 1;
    }
  }

  Search S;
  // ---- CNF ----
  {
    std::ifstream in(argv[1]);
    std::string tok;
    int nclauses = 0;
    std::vector<int> cur;
    while (in >> tok) {
      if (tok == "c") { std::string rest; std::getline(in, rest); continue; }
      if (tok == "p") { std::string cnf; in >> cnf >> S.nvars >> nclauses; continue; }
      int lit = atoi(tok.c_str());
      if (lit == 0) { S.clauses.push_back(cur); cur.clear(); }
      else cur.push_back(lit);
    }
  }
  std::vector<int> watched, block;
  {
    std::ifstream in(argv[2]);
    std::string line;
    if (std::getline(in, line)) watched = parse_ints(line);
    if (std::getline(in, line)) block = parse_ints(line);
  }
  std::vector<int> rule_bits;
  {
    std::ifstream in(argv[3]);
    std::string line;
    std::getline(in, line);
    for (char ch : line)
      if (ch == '0' || ch == '1') rule_bits.push_back(ch - '0');
    if (rule_bits.size() != 512) {
      emit("INFO bad rule table");
      return 1;
    }
    std::getline(in, line);
    S.order = parse_ints(line);
    while (std::getline(in, line)) {
      std::istringstream ls(line);
      std::array<int, NPOS> t;
      std::string w;
      int n = 0;
      while (n < NPOS && ls >> w) {
        if (w == "T") t[n] = TRUE_LIT;
        else if (w == "F") t[n] = FALSE_LIT;
        else {
          t[n] = atoi(w.c_str());
          if (abs(t[n]) > S.nvars) S.nvars = abs(t[n]);
        }
        n++;
      }
      if (n == NPOS) S.trans.push_back(t);
    }
  }
  for (int v : watched) S.nvars = std::max(S.nvars, v);
  for (int l : S.order) S.nvars = std::max(S.nvars, abs(l));

  S.build_table(rule_bits);
  int nv = S.nvars;
  S.value.assign(nv + 1, 0);
  S.level.assign(nv + 1, 0);
  S.watches.assign(2 * (nv + 1) + 2, {});
  S.occurs.assign(nv + 1, {});
  S.is_block.assign(nv + 1, 0);
  for (int v : block)
    if (v > 0 && v <= nv) S.is_block[v] = 1;
  for (size_t id = 0; id < S.trans.size(); id++) {
    for (int i = 0; i < NPOS; i++) {
      int l = S.trans[id][i];
      if (l == TRUE_LIT || l == FALSE_LIT) continue;
      auto &o = S.occurs[abs(l)];
      if (o.empty() || o.back() != (int) id) o.push_back((int) id);
    }
  }
  size_t given = 0; // decision-order entries that are cells (the rest are auxiliaries)
  // decision order: given literals, then every remaining variable (false first)
  {
    std::vector<char> seen(nv + 1, 0);
    std::vector<int> order;
    for (int l : S.order)
      if (l != 0 && !seen[abs(l)]) { seen[abs(l)] = 1; order.push_back(l); }
    given = order.size();
    for (int v = 1; v <= nv; v++)
      if (!seen[v]) order.push_back(-v);
    S.order = order;
  }

  char info[160];
  snprintf(info, sizeof info, "INFO backtracker: %d variables, %zu clauses, %zu transitions, watching %zu",
           nv, S.clauses.size(), S.trans.size(), watched.size());
  emit(info);

  std::thread reader([] {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line == "pause") g_pause = true;
      else if (line == "resume") g_pause = false;
      else if (line == "stop") { g_stop = true; g_pause = false; }
      else if (line == "snap") g_snap = true;
    }
    g_stop = true;
    g_pause = false;
  });
  reader.detach();

  std::atomic<long long> total_decisions{0}, total_backtracks{0}, found{0};
  auto snapshot = [&](const Search &W) {
    std::string s;
    s.reserve(watched.size() + 64);
    char head[96];
    snprintf(head, sizeof head, "SNAP %.3f %lld %lld ", elapsed(), total_backtracks.load(), total_decisions.load());
    s += head;
    for (int v : watched) {
      int x = (v > 0 && v <= nv) ? W.value[v] : 0;
      bool root = x != 0 && W.level[v] == 0;
      s += x > 0 ? (root ? 't' : '1') : x < 0 ? (root ? 'f' : '0') : '.';
    }
    emit(s);
  };
  auto finish = [&](const char *what) {
    char buf[96];
    if (what) {
      snprintf(buf, sizeof buf, "%s %.3f", what, elapsed());
      emit(buf);
    }
    snprintf(buf, sizeof buf, "DONE %.3f %lld", elapsed(), found.load());
    emit(buf);
    fflush(stdout);
    _exit(0);
  };

  // ---- level 0: unit and empty clauses, constant-only transitions, watches ----
  bool ok = true;
  for (size_t cid = 0; cid < S.clauses.size() && ok; cid++) {
    auto &c = S.clauses[cid];
    if (c.empty()) ok = false;
    else if (c.size() == 1) ok = S.assign(c[0]);
    else {
      S.watches[Search::lit_index(c[0])].push_back((int) cid);
      S.watches[Search::lit_index(c[1])].push_back((int) cid);
    }
  }
  for (size_t id = 0; id < S.trans.size() && ok; id++) ok = S.check_transition((int) id);
  // Clauses whose watched literals were made false before watching began need a pass:
  // re-propagate every assigned literal from the start.
  S.qhead = 0;
  if (ok) ok = S.propagate();
  if (ok) {
    // binary/longer clauses may have been unit or falsified by root assignments that
    // were made before they were watched; check them all once
    for (size_t cid = 0; cid < S.clauses.size() && ok; cid++) {
      auto &c = S.clauses[cid];
      if (c.size() < 2) continue;
      int unknown = 0, last = 0;
      bool sat = false;
      for (int l : c) {
        int v = S.lit_value(l);
        if (v == 1) { sat = true; break; }
        if (v == 0) { unknown++; last = l; }
      }
      if (sat) continue;
      if (unknown == 0) ok = false;
      else if (unknown == 1) { ok = S.assign(last) && S.propagate(); }
    }
  }
  if (!ok) finish("UNSAT");

  if (probe) {
    // Failed-literal probing at the root, repeated while it keeps fixing cells
    long long fixed = 0;
    for (int round = 0; round < 3; round++) {
      long long before = fixed;
      for (size_t k = 0; k < given && !g_stop; k++) {
        int lit = S.order[k];
        if (S.value[abs(lit)] != 0) continue;
        bool works[2];
        for (int side = 0; side < 2; side++) {
          size_t mark = S.trail.size();
          S.trail_lim.push_back(mark);
          S.flipped.push_back(0);
          S.assign(side ? -lit : lit);
          works[side] = S.propagate();
          S.undo_to(mark);
          S.trail_lim.pop_back();
          S.flipped.pop_back();
        }
        if (!works[0] && !works[1]) finish("UNSAT");
        if (!works[0] || !works[1]) {
          fixed++;
          if (!S.assign(works[0] ? lit : -lit) || !S.propagate()) finish("UNSAT");
        }
      }
      if (fixed == before) break;
    }
    snprintf(info, sizeof info, "INFO probing fixed %lld cells", fixed);
    emit(info);
  }

  // ---- search: the tree is split at depth `depth` into 2^depth parts, which the
  // threads take from a shared counter. In a part, the first `depth` decisions follow
  // the part number's bits and are never flipped.
  int depth = 0;
  if (threads > 1)
    while ((1 << depth) < threads * 16 && depth < 20) depth++;
  const long long nparts = 1LL << depth;
  std::atomic<long long> next_part{0};
  std::atomic<int> viewer{-1};  // the thread whose state the live view shows
  std::atomic<long long> last_snap_ms{0};
  std::atomic<bool> pause_announced{false};

  auto worker = [&](int tid) {
    unsigned counter = 0;
    long long part;
    while ((part = next_part.fetch_add(1)) < nparts) {
      Search W = S; // private copy of the root state
      int expected = -1;
      viewer.compare_exchange_strong(expected, tid);
      long long reported_d = 0, reported_b = 0;
      while (true) {
        if ((++counter & 1023) == 0 || g_snap || g_pause || g_stop) {
          total_decisions += W.decisions - reported_d;
          total_backtracks += W.backtracks - reported_b;
          reported_d = W.decisions;
          reported_b = W.backtracks;
          if (g_stop) finish("STOPPED");
          if (g_pause) {
            if (!pause_announced.exchange(true)) {
              snapshot(W);
              emit("INFO paused");
            }
            while (g_pause && !g_stop) {
              if (viewer == tid && g_snap.exchange(false)) snapshot(W);
              std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
            if (g_stop) finish("STOPPED");
            if (pause_announced.exchange(false)) emit("INFO resumed");
          }
          if (viewer == tid) {
            long long now = (long long) (elapsed() * 1000);
            if (g_snap.exchange(false) || now - last_snap_ms >= interval_ms) {
              last_snap_ms = now;
              snapshot(W);
            }
          }
        }

        int lit = W.next_decision();
        if (lit == 0) {
          // every variable assigned without conflict: a solution. If it was reached in
          // fewer than `depth` decisions, several parts contain it; only the one whose
          // remaining bits are zero reports it.
          int d = (int) W.trail_lim.size();
          bool mine = d >= depth || (part & ((1LL << (depth - d)) - 1)) == 0;
          if (mine && (max_solutions == 0 || found < max_solutions)) {
            long long n = ++found;
            if (max_solutions != 0 && n > max_solutions) break; // another thread got there first
            std::string out;
            char head[64];
            snprintf(head, sizeof head, "SOL %lld %.3f ", n, elapsed());
            out += head;
            for (int v : watched) out += (v > 0 && v <= nv && W.value[v] > 0) ? '1' : '0';
            emit(out);
            if (max_solutions != 0 && n >= max_solutions) finish(nullptr);
          }
          W.skip_equivalent_solutions();
          if (!W.backtrack()) break;
          continue;
        }
        int d = (int) W.trail_lim.size();
        bool fixed_decision = false;
        if (d < depth) {
          if (W.is_block[abs(lit)]) {
            // a fixed decision of this part: follow its bit, no alternative
            if ((part >> (depth - 1 - d)) & 1) lit = -lit;
            fixed_decision = true;
          } else if (part & ((1LL << (depth - d)) - 1)) {
            // Only auxiliary variables are left, so the pattern is already decided and
            // every part sharing the bits so far holds the same solutions. The part whose
            // remaining bits are zero searches this subtree; the others skip it.
            break;
          }
        }
        W.decisions++;
        W.trail_lim.push_back(W.trail.size());
        W.flipped.push_back(fixed_decision ? 1 : 0);
        W.assign(lit);
        if (!W.propagate()) {
          if (!W.backtrack()) break;
        }
      }
      total_decisions += W.decisions - reported_d;
      total_backtracks += W.backtracks - reported_b;
      expected = tid;
      viewer.compare_exchange_strong(expected, -1);
    }
  };

  if (threads <= 1) worker(0);
  else {
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, t);
    for (auto &t : pool) t.join();
  }
  finish("UNSAT");
}
