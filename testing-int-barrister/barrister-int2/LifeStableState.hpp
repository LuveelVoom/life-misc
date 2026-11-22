#pragma once

// Generic still-life solver for the INT port of Barrister.
//
// The stable background is a three-valued grid (state, unknown). Its only
// constraint is the rule itself: at every cell, the rule applied to the 3x3
// neighbourhood gives back the centre. gIntRule.still (IntRule.hpp) holds, for
// every on/off/unknown 3x3, whether some completion satisfies that cell's
// constraint and which unknown cells are forced. Propagation applies it until
// nothing changes. This replaces upstream Barrister's Life-specific neighbour
// count options (live2/live3/dead0..dead6), so it is exact for any INT or MAP
// rule without B0.
//
// Every window containing a known ON cell is checked, whether or not it still
// has unknowns. Any window without one is satisfied by turning its unknowns
// off (no B0), so no inconsistency can be missed. Forced cells in windows without a known ON cell (possible in rules
// with B1/B2) are not deduced; that only weakens pruning.

#include <chrono>
#include <limits>

#include "LifeAPI.h"
#include "Parsing.hpp"
#include "LifeHistoryState.hpp"

struct PropagateResult {
  bool consistent;
  bool changed;
};

enum struct CompletionResult {
  COMPLETED,
  INCONSISTENT,
  TIMEOUT,
};

class LifeStableState {
public:
  LifeState state;
  LifeState unknown;
  // Cells that are ON or next to a known ON cell: the catalyst's zone of
  // influence. Name kept from upstream, where it meant "dead with 0 live
  // neighbours ruled out"; Barrister uses it as that zone.
  LifeState dead0;
  // Cells set since the last successful Propagate; only windows around them
  // are re-checked. fullCheck forces a check of every window (after state or
  // unknown were assigned directly).
  LifeState dirty;
  bool fullCheck = true;

  bool operator==(const LifeStableState &o) const { return state == o.state && unknown == o.unknown; }

  void UpdateZOI() { dead0 = state.ZOI(); }

  void Move(int x, int y) {
    state.Move(x, y);
    unknown.Move(x, y);
    fullCheck = true;
    UpdateZOI();
  }
  void Move(std::pair<int, int> vec) { Move(vec.first, vec.second); }

  void SetOn(const LifeState &which) { state |= which; unknown &= ~which; dirty |= which; UpdateZOI(); }
  void SetOff(const LifeState &which) { state &= ~which; unknown &= ~which; dirty |= which; UpdateZOI(); }
  void SetOn(std::pair<int, int> cell) { state.Set(cell); unknown.Erase(cell); dirty.Set(cell); UpdateZOI(); }
  void SetOff(std::pair<int, int> cell) { state.Erase(cell); unknown.Erase(cell); dirty.Set(cell); UpdateZOI(); }

  // Upstream kept count options in sync with state/unknown; nothing to do here.
  PropagateResult SynchroniseStateKnown() { return {true, false}; }
  PropagateResult SynchroniseStateKnown(std::pair<int, int>) { return {true, false}; }

  PropagateResult Propagate();
  PropagateResult PropagateStrip(unsigned) { return Propagate(); }

  // Value of each unknown cell is tried both ways; a value that is
  // inconsistent forces the other.
  LifeState Vulnerable() const { return unknown & state.ZOI(); }
  PropagateResult TestUnknown(std::pair<int, int> cell);
  PropagateResult TestUnknowns(const LifeState &cells);

  LifeStableState Join(const LifeStableState &other) const;
  LifeStableState Graft(const LifeStableState &other) const;
  LifeStableState ClearUnmodified() const;
  LifeState Differences(const LifeStableState &other) const {
    return (state ^ other.state) | (unknown ^ other.unknown);
  }

  CompletionResult CompleteStableStep(std::chrono::system_clock::time_point &timeLimit, bool minimise,
                                      bool useSeed, const LifeState &seed, unsigned &maxPop, LifeState &best);
  std::pair<CompletionResult, LifeState> CompleteStable(unsigned timeout, bool minimise);

  bool CompatibleWith(const LifeState &desired) const {
    return (~unknown & (state ^ desired)).IsEmpty();
  }
  bool CompatibleWith(const LifeStableState &desired) const {
    return (~unknown & ~desired.unknown & (state ^ desired.state)).IsEmpty();
  }

  void SanityCheck() const {
#ifdef DEBUG
    assert((state & unknown).IsEmpty());
#endif
  }

  std::string RLE() const { return LifeBellmanRLEFor(state, unknown | state); }
  std::string RLEWHeader() const { return "x = 0, y = 0, rule = LifeBellman\n" + RLE(); }
  friend std::ostream &operator<<(std::ostream &os, LifeStableState const &self) {
    return os << self.RLEWHeader();
  }
};

// Neighbourhood positions in IntIndex9 bit order: bit 0 SE, 1 S, 2 SW, 3 E,
// 4 W, 5 NE, 6 N, 7 NW, 8 centre (x = column, y = bit, y down).
static constexpr int kPosDX[9] = {1, 0, -1, 1, -1, 1, 0, -1, 0};
static constexpr int kPosDY[9] = {1, 1, 1, 0, 0, -1, -1, -1, 0};

inline unsigned StableIndex9(const LifeState &s, unsigned x, unsigned y) {
  return IntIndex9(s[(x + N - 1) % N], s[x], s[(x + 1) % N], y);
}

PropagateResult LifeStableState::Propagate() {
  bool anyChange = false;
  LifeState todo = fullCheck ? state.ZOI() : (dirty.ZOI() & state.ZOI());
  while (!todo.IsEmpty()) {
    LifeState changedCells;
    for (unsigned x = 0; x < N; x++) {
      uint64_t col = todo[x];
      while (col) {
        unsigned y = __builtin_ctzll(col);
        col &= col - 1;
        unsigned on9 = StableIndex9(state, x, y);
        unsigned unk9 = StableIndex9(unknown, x, y);
        uint32_t v = gIntRule.still[(unk9 << 9) | on9];
        if (!(v >> 31)) {
          fullCheck = true; // leave the state to be fully re-checked
          return {false, true};
        }
        uint32_t fon = v & 511u, foff = (v >> 9) & 511u;
        for (unsigned b = 0; b < 9; b++) {
          if (!((fon | foff) >> b & 1)) continue;
          unsigned cx = (x + N + kPosDX[b]) % N, cy = (y + 64 + kPosDY[b]) & 63;
          if (fon >> b & 1) state.Set(cx, cy);
          unknown.Erase(cx, cy);
          changedCells.Set(cx, cy);
        }
      }
    }
    if (changedCells.IsEmpty())
      break;
    anyChange = true;
    todo = changedCells.ZOI() & state.ZOI();
  }
  dirty = LifeState();
  fullCheck = false;
  UpdateZOI();
  return {true, anyChange};
}

// Join of two propagated, consistent refinements of the same state: keeps the
// cells both agree on. Those are implied facts, so the result is consistent;
// only the cells that became known need their windows re-checked.
LifeStableState LifeStableState::Join(const LifeStableState &other) const {
  LifeStableState result;
  result.unknown = unknown | other.unknown | (state ^ other.state);
  result.state = state & other.state & ~result.unknown;
  result.dirty = ~result.unknown & (unknown | other.unknown);
  result.fullCheck = fullCheck || other.fullCheck;
  result.UpdateZOI();
  return result;
}

// Cells unknown here but known in other take other's value.
LifeStableState LifeStableState::Graft(const LifeStableState &other) const {
  LifeStableState result = *this;
  LifeState take = unknown & ~other.unknown;
  result.state = (state & ~take) | (other.state & take);
  result.unknown = unknown & ~take;
  result.fullCheck = true;
  result.UpdateZOI();
  return result;
}

LifeStableState LifeStableState::ClearUnmodified() const {
  LifeStableState result = *this;
  LifeState toClear = unknown & ~(dead0.ZOI());
  result.unknown &= ~toClear;
  result.state &= ~toClear;
  result.fullCheck = true;
  result.UpdateZOI();
  return result;
}

PropagateResult LifeStableState::TestUnknown(std::pair<int, int> cell) {
  LifeStableState onSearch = *this;
  onSearch.SetOn(cell);
  auto onResult = onSearch.Propagate();

  LifeStableState offSearch = *this;
  offSearch.SetOff(cell);
  auto offResult = offSearch.Propagate();

  if (!onResult.consistent && !offResult.consistent)
    return {false, false};
  if (!onResult.consistent) {
    *this = offSearch;
    return {true, true};
  }
  if (!offResult.consistent) {
    *this = onSearch;
    return {true, true};
  }
  if (onResult.changed && offResult.changed) {
    LifeStableState joined = onSearch.Join(offSearch);
    joined.dirty = unknown & ~joined.unknown;   // newly known relative to this state
    joined.fullCheck = false;   // both branches were fully checked and consistent
    bool change = joined != *this;
    *this = joined;
    return {true, change};
  }
  return {true, false};
}

PropagateResult LifeStableState::TestUnknowns(const LifeState &cells) {
  LifeState remainingCells = cells & unknown;
  bool anyChanges = false;
  while (!remainingCells.IsEmpty()) {
    auto cell = remainingCells.FirstOn();
    remainingCells.Erase(cell);
    auto result = TestUnknown(cell);
    if (!result.consistent)
      return {false, false};
    anyChanges = anyChanges || result.changed;
    remainingCells &= unknown;
  }
  return {true, anyChanges};
}

// Branches on unknown cells next to a known ON cell (OFF first). When none are
// left, all remaining unknowns are tried OFF and checked exactly; if that
// fails, cells at distance 2 are branched on too. Unknowns at distance >= 3
// from every ON cell can always be OFF (no B0).
CompletionResult LifeStableState::CompleteStableStep(
    std::chrono::system_clock::time_point &timeLimit, bool minimise,
    bool useSeed, const LifeState &seed, unsigned &maxPop, LifeState &best) {
  if (std::chrono::system_clock::now() > timeLimit)
    return CompletionResult::TIMEOUT;

  if (!Propagate().consistent)
    return CompletionResult::INCONSISTENT;

  if (state.GetPop() >= maxPop)
    return CompletionResult::COMPLETED;

  LifeState settable = unknown & state.ZOI();
  if (settable.IsEmpty()) {
    // Every unknown cell is at distance >= 2 from all ON cells. Turning them
    // all OFF usually works; check it exactly, else branch on distance 2.
    LifeStableState allOff = *this;
    allOff.state &= ~allOff.unknown;
    allOff.unknown = LifeState();
    allOff.fullCheck = true;
    if (allOff.Propagate().consistent) {
      best = state;
      maxPop = state.GetPop();
      return CompletionResult::COMPLETED;
    }
    settable = unknown & state.ZOI().ZOI();
    if (settable.IsEmpty())
      return CompletionResult::INCONSISTENT;
  }

  if (useSeed) {
    LifeState seedZOI = seed;
    while ((settable & seedZOI).IsEmpty())
      seedZOI = seedZOI.ZOI();
    settable &= seedZOI;
  }

  std::pair<int, int> cell = (Vulnerable() & settable).FirstOn();
  if (cell.first == -1)
    cell = settable.FirstOn();

  {
    LifeStableState next = *this;
    next.SetOff(cell);
    auto result = next.CompleteStableStep(timeLimit, minimise, useSeed, seed, maxPop, best);
    if (result == CompletionResult::TIMEOUT)
      return CompletionResult::TIMEOUT;
    if (!minimise && result == CompletionResult::COMPLETED)
      return CompletionResult::COMPLETED;
  }
  SetOn(cell);
  return CompleteStableStep(timeLimit, minimise, useSeed, seed, maxPop, best);
}

std::pair<CompletionResult, LifeState> LifeStableState::CompleteStable(unsigned timeout, bool minimise) {
  if (state.IsEmpty())
    return {CompletionResult::COMPLETED, LifeState()};
  if (unknown.IsEmpty())
    return {CompletionResult::COMPLETED, state};

  LifeState best;
  unsigned maxPop = std::numeric_limits<int>::max();
  auto timeLimit = std::chrono::system_clock::now() + std::chrono::seconds(timeout);
  CompletionResult result = CompletionResult::INCONSISTENT;

  // First look for a completion with a small bounding box
  LifeState searchArea = state;
  while (!(unknown & ~searchArea).IsEmpty()) {
    searchArea = searchArea.ZOI();
    LifeStableState copy = *this;
    copy.unknown &= searchArea;
    copy.fullCheck = true;
    copy.UpdateZOI();
    result = copy.CompleteStableStep(timeLimit, minimise, false, state, maxPop, best);
    if (best.GetPop() > 0 || std::chrono::system_clock::now() > timeLimit)
      break;
  }

  if (best.IsEmpty())
    return {result == CompletionResult::TIMEOUT ? CompletionResult::TIMEOUT : CompletionResult::INCONSISTENT,
            LifeState()};

  if (minimise) {
    LifeStableState copy = *this;
    copy.unknown &= searchArea.BigZOI();
    copy.fullCheck = true;
    copy.CompleteStableStep(timeLimit, minimise, true, state | best, maxPop, best);
  }
  return {CompletionResult::COMPLETED, best};
}
