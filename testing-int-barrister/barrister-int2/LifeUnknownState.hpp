#pragma once

#include <bit>

#include "LifeAPI.h"
#include "Bits.hpp"
#include "LifeStableState.hpp"
#include "Transition.hpp"
#include "LifeHistoryState.hpp"

class LifeUnknownState {
public:
  LifeState state;
  LifeState unknown;
  LifeState unknownStable; // Equal to the stable state

  bool operator==(const LifeUnknownState &b) const {
    return state == b.state && unknown == b.unknown && unknownStable == b.unknownStable;
  }
  bool operator!=(const LifeUnknownState &b) const {
    return !(*this == b);
  }

  LifeUnknownState StepMaintaining(const LifeStableState &stable) const;
  std::tuple<uint64_t, uint64_t, uint64_t> StepMaintainingColumn(const LifeStableState &stable, int i) const;
  std::tuple<std::array<uint64_t, 4>, std::array<uint64_t, 4>, std::array<uint64_t, 4>> StepMaintainingStrip(const LifeStableState &stable, int i) const;

  // bool CanCleanlyAdvance(const LifeStableState &stable) const;
  LifeState ActiveComparedTo(const LifeStableState &stable) const;
  LifeState ChangesComparedTo(const LifeUnknownState &prev) const;
  void TransferStable(const LifeStableState &stable);
  void TransferStable(const LifeStableState &stable, std::pair<int, int> cell);

  void SetKnown(std::pair<int, int> cell, bool value, bool stable) {
    if (stable) {
      unknown.Set(cell);
      unknownStable.Set(cell);
    } else {
      state.Set(cell, value);
      unknown.Erase(cell);
      unknownStable.Erase(cell);
    }
  }

  void SetTransitionResult(std::pair<int, int> cell, Transition transition) {
    switch (transition) {
    case Transition::OFF_TO_OFF:
    case Transition::ON_TO_OFF:
      SetKnown(cell, false, false);
      break;
    case Transition::OFF_TO_ON:
    case Transition::ON_TO_ON:
      SetKnown(cell, true, false);
      break;
    case Transition::STABLE_TO_STABLE:
      SetKnown(cell, false, true);
      break;
    default:
      break;
    }
  }

  void SetTransitionPrev(std::pair<int, int> cell, Transition transition) {
    switch (transition) {
    case Transition::OFF_TO_OFF:
    case Transition::OFF_TO_ON:
      SetKnown(cell, false, false);
      break;
    case Transition::ON_TO_OFF:
    case Transition::ON_TO_ON:
      SetKnown(cell, true, false);
      break;
    case Transition::STABLE_TO_STABLE:
      SetKnown(cell, false, true);
      break;
    default:
      break;
    }
  }

  Transition UnperturbedTransitionFor(std::pair<int, int> cell) const {
    bool next = state.StepFor(cell);
    switch (state.Get(cell) << 1 | next) {
    case 0b00: return Transition::OFF_TO_OFF;
    case 0b01: return Transition::OFF_TO_ON;
    case 0b10: return Transition::ON_TO_OFF;
    case 0b11: return Transition::ON_TO_ON;
    default: return Transition::IMPOSSIBLE;
    }
  }

  // Would this transition be an interaction with the active pattern?
  bool TransitionIsPerturbation(std::pair<int, int> cell, Transition transition) const {
    return transition != Transition::STABLE_TO_STABLE && transition != UnperturbedTransitionFor(cell);
  }

  // Would this transition be a active cell?
  bool TransitionIsActive(std::pair<int, int> cell, Transition transition) const {
    return (transition == Transition::ON_TO_OFF ||
            transition == Transition::OFF_TO_ON) &&
           TransitionIsPerturbation(cell, transition);
  }

  std::tuple<bool, bool, bool> StepMaintainingFor(const LifeStableState &stable,
                                                  std::pair<int, int> cell) const {
    auto [next, nextUnknown, nextUnknownStable] = StepMaintainingColumn(stable, cell.first);

    int y = cell.second;
    bool cellNext    = (next                    & (1ULL << y)) >> y;
    bool cellUnknown = (nextUnknown             & (1ULL << y)) >> y;
    bool cellUnknownStable = (nextUnknownStable & (1ULL << y)) >> y;
    return {cellNext, cellUnknown, cellUnknownStable};
  }

  Transition TransitionsFor(const LifeStableState &stable, std::pair<int, int> cell) const {
    bool prevUnknown = unknown.Get(cell) && !(unknownStable.Get(cell) && !stable.unknown.Get(cell));
    bool prevState = state.Get(cell) || (unknownStable.Get(cell) && stable.state.Get(cell));
    auto [nextState, nextUnknown, nextStableUnknown] = StepMaintainingFor(stable, cell);

    Transition transitions = Transition::ANY;
    if (!prevUnknown) {
      if (prevState)
        transitions &= Transition::ON_TO_OFF | Transition::ON_TO_ON;
      else
        transitions &= Transition::OFF_TO_OFF | Transition::OFF_TO_ON;
    }
    if (!nextUnknown) {
      if (nextState)
        transitions &= Transition::OFF_TO_ON | Transition::ON_TO_ON;
      else
        transitions &= Transition::OFF_TO_OFF | Transition::ON_TO_OFF;
    }
    return transitions;
  }
  
  LifeHistoryState ToHistory() const {
    return LifeHistoryState(state, unknown, unknownStable, LifeState());
  }

  friend std::ostream& operator<<(std::ostream& os, LifeUnknownState const& self) {
    return os << self.ToHistory().RLEWHeader();
  }

  void SanityCheck(const LifeStableState &stable) {
#ifdef DEBUG
    assert((unknown & state).IsEmpty());
    assert((~unknown & unknownStable).IsEmpty());
    assert((~unknown & stable.unknown).IsEmpty());
    // assert((unknownStable & ~stable.unknown).IsEmpty());
#endif
  }
};

// INT port: StepMaintaining is rewritten per cell instead of the Life
// bitsliced circuits. For each cell:
//  1. no unknown cell in its 3x3: exact rule step;
//  2. otherwise, if its 3x3 matches the stable background exactly (every cell
//     known-equal to a known stable cell, or a stable unknown): it keeps the
//     stable value (unknown-stable if the stable cell is unknown). This relies
//     on the background being a still life of the rule, which the final
//     exact check verifies;
//  3. otherwise, an unknown centre becomes an active unknown;
//  4. otherwise a three-valued lookup over all completions of the unknowns.
// Unlike the original, this ignores the stable options of the cells (Life
// count reasoning), so it is weaker but does not assume Life's rule.
static inline uint64_t IntDiffColumn(const LifeUnknownState &s, const LifeStableState &stable, unsigned c) {
  return ((s.state[c] ^ stable.state[c]) & ~s.unknown[c]) |
         (s.unknown[c] & ~s.unknownStable[c]) |
         (~s.unknown[c] & stable.unknown[c]);
}

static inline std::tuple<uint64_t, uint64_t, uint64_t>
IntStepMaintainingColumn(const LifeUnknownState &s, const LifeStableState &stable, int i) {
  unsigned l = (i + N - 1) % N, r = (i + 1) % N;
  uint64_t L = s.state[l], C = s.state[i], R = s.state[r];
  uint64_t UL = s.unknown[l], UC = s.unknown[i], UR = s.unknown[r];

  uint64_t ucol = UL | UC | UR;
  uint64_t unkZ = std::rotl(ucol, 1) | ucol | std::rotr(ucol, 1);

  uint64_t col = L | C | R;
  uint64_t liveZ = std::rotl(col, 1) | col | std::rotr(col, 1);

  // The stable background is a still life of the rule (generic solver), so a
  // cell whose whole 3x3 equals the background keeps its background value.
  uint64_t dcol = IntDiffColumn(s, stable, l) | IntDiffColumn(s, stable, i) | IntDiffColumn(s, stable, r);
  uint64_t diffZ = std::rotl(dcol, 1) | dcol | std::rotr(dcol, 1);
  uint64_t restore = ~diffZ;

  uint64_t next_on = (restore & stable.state[i] & ~stable.unknown[i]) | IntStepBits(L, C, R, diffZ & ~unkZ & liveZ);
  uint64_t next_unknown = restore & stable.unknown[i];
  uint64_t next_unknown_stable = restore & stable.unknown[i];

  if (unkZ & diffZ) {

    next_unknown |= unkZ & diffZ & UC;

    uint64_t todo = unkZ & diffZ & ~UC;
    while (todo) {
      unsigned y = __builtin_ctzll(todo);
      todo &= todo - 1;
      unsigned on9 = IntIndex9(L, C, R, y);
      unsigned unk9 = IntIndex9(UL, UC, UR, y);
      uint8_t t = gIntRule.tri[(unk9 << 9) | on9];
      if (t == 1) next_on |= 1ULL << y;
      else if (t == 2) next_unknown |= 1ULL << y;
    }
  }
  return {next_on, next_unknown, next_unknown_stable};
}

LifeUnknownState LifeUnknownState::StepMaintaining(const LifeStableState &stable) const {
  LifeUnknownState result{LifeState(false), LifeState(false), LifeState(false)};
  for (int i = 0; i < N; i++) {
    auto [on, unk, unkStable] = IntStepMaintainingColumn(*this, stable, i);
    result.state[i] = on;
    result.unknown[i] = unk;
    result.unknownStable[i] = unkStable;
  }
  return result;
}

std::tuple<uint64_t, uint64_t, uint64_t> LifeUnknownState::StepMaintainingColumn(const LifeStableState &stable, int i) const {
  return IntStepMaintainingColumn(*this, stable, i);
}

std::tuple<std::array<uint64_t, 4>,
           std::array<uint64_t, 4>,
           std::array<uint64_t, 4>>
LifeUnknownState::StepMaintainingStrip(const LifeStableState &stable, int column) const {
  std::array<uint64_t, 4> resultState;
  std::array<uint64_t, 4> resultUnknown;
  std::array<uint64_t, 4> resultUnknownStable;
  // Same columns as GetStrip<4>(column): column-1 .. column+2
  for (int k = 0; k < 4; k++) {
    int c = (column + k - 1 + N) % N;
    std::tie(resultState[k], resultUnknown[k], resultUnknownStable[k]) =
        IntStepMaintainingColumn(*this, stable, c);
  }
  return {resultState, resultUnknown, resultUnknownStable};
}

LifeState LifeUnknownState::ActiveComparedTo(const LifeStableState &stable) const {
  return ~unknown & ~stable.unknown & (stable.state ^ state);
}

void LifeUnknownState::TransferStable(const LifeStableState &stable) {
  LifeState updated = unknownStable & ~stable.unknown;
  state |= stable.state & updated;
  unknown &= ~updated;
  unknownStable &= ~updated;
}

void LifeUnknownState::TransferStable(const LifeStableState &stable, std::pair<int, int> cell) {
  bool updated = unknownStable.Get(cell) && !stable.unknown.Get(cell);
  if (updated) {
    if(stable.state.Get(cell))
      state.Set(cell);
    unknown.Erase(cell);
    unknownStable.Erase(cell);
  }
}

LifeState
LifeUnknownState::ChangesComparedTo(const LifeUnknownState &prev) const {
  return (state ^ prev.state) & ~unknown & ~prev.unknown;
}
