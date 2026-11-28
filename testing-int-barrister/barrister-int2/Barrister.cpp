#include <cassert>
#include <cmath>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>

#include "toml/toml.hpp"

#include "Bits.hpp"
#include "LifeAPI.h"
#include "LifeStableState.hpp"
#include "LifeUnknownState.hpp"
#include "RotorDescription.hpp"
#include "Params.hpp"
#include "Parsing.hpp"

// Idea:
//
// * Calculate the frontier

// * Go through it a single cell at a time

// * For each cell, determine which of the cases (ON to OFF, OFF to
// ON, STABLE to STABLE) are allowed. If there is only one then set it
// and remove the cell from the frontier

// * If we set any cell of the frontier, then propagate stable and
// start over

// * Once there are only true choices left, branch on a cell in the
// frontier (earliest first?)

// After branching, we will likely get some new frontier cells, and
// some previous frontier cells might become settable. Maybe we should
// check the existing frontier cells for any settable ones quickly,
// before doing a full recalculation. (Maybe only if we have made the
// new cell active instead of stable)

// It might also be useful to calculate the "semi-frontier" while
// calculating the frontier, and using that information somehow when
// choosing which cell to branch on. Perhaps something like: branch on
// the frontier cell whose ZOI (over all generations) touches the most
// semi-frontier cells

// The UnknownStep could be made more intelligent by handling unknown
// active cells in the neighbourhood. E.g., if we are a DEAD6 cell
// then we will stay dead.

// It may not be worth the cost: if there is an unknown active cell
// then the current cell is likely to be swamped in the next
// generation or two

const unsigned maxFrontierGens = 6;
const unsigned maxBranchFastCount = 1;
const unsigned maxCalculateRounds = 1;

const unsigned maxCellActiveWindowGens = 0;
const unsigned maxCellActiveStreakGens = 0;

struct Solution {
  LifeState state;
  LifeState completed;
  LifeStableState stable;
  LifeStableState interactionStable;
  LifeState stator;

  unsigned interactionGen;
  unsigned recoveryGen;

  CompletionResult completionResult;

  bool operator==(const Solution&) const = default; // I don't really know why I need to say this

  // Doesn't have to be fast
  auto operator<=>(const Solution &other) const {
    if (auto c = interactionGen <=> other.interactionGen; c != 0) return c;
    if (auto c = stable.state.GetPop() <=> other.stable.state.GetPop(); c != 0) return c;
    if (auto c = recoveryGen <=> other.recoveryGen; c != 0) return c;
    return stable.state.GetHash() <=> other.stable.state.GetHash();
  }
};

Transition AllowedTransitions(bool state, bool unknownstable, bool stablestate,
                              bool forcedInactive, bool forcedUnchanging, bool inzoi, Transition unperturbed) {
  auto result = Transition::ANY & ~Transition::STABLE_TO_STABLE;

  // If current state is known, remove options with the wrong previous state
  if (!unknownstable && state)
    result &= Transition::ON_TO_OFF | Transition::ON_TO_ON;
  if (!unknownstable && !state)
    result &= Transition::OFF_TO_OFF | Transition::OFF_TO_ON;

  if (forcedInactive && inzoi) {
    if (unknownstable)
      result &= ~(Transition::OFF_TO_ON | Transition::ON_TO_OFF);
    if (!unknownstable && stablestate)
      result &= ~(Transition::OFF_TO_OFF | Transition::ON_TO_OFF);
    if (!unknownstable && !stablestate)
      result &= ~(Transition::OFF_TO_ON | Transition::ON_TO_ON);
  }

  if (forcedInactive && !inzoi) {
    result &= unperturbed | Transition::OFF_TO_OFF | Transition::ON_TO_ON;
  }

  if (forcedUnchanging && inzoi)
    result &= Transition::OFF_TO_OFF | Transition::ON_TO_ON;

  return result;
}

struct FrontierGeneration {
  LifeUnknownState state;
  LifeUnknownState next;
  LifeState frontierCells;
  LifeState active;
  LifeState changes;
  LifeState forcedInactive;
  LifeState forcedUnchanging;
  unsigned gen;

  void SetTransition(std::pair<int, int> cell, Transition transition) {
    state.SetTransitionPrev(cell, transition);
    next.SetTransitionResult(cell, transition);
  }

  std::string RLE() const {
    LifeHistoryState history(next.state, next.unknown & ~next.unknownStable,
                             next.unknownStable, LifeState());
    return history.RLEWHeader();
  }

  Transition AllowedTransitions(const LifeStableState &stable,
                              std::pair<int, int> cell) const {
    auto allowedTransitions = ::AllowedTransitions(
        state.state.Get(cell), stable.unknown.Get(cell), stable.state.Get(cell),
        forcedInactive.Get(cell), forcedUnchanging.Get(cell),
        stable.dead0.Get(cell), state.UnperturbedTransitionFor(cell));

    return allowedTransitions;
  }
};

class SearchState {
public:
  LifeStableState stable;
  FrontierGeneration frontier;

  LifeState everActive;

  LifeCountdown<maxCellActiveWindowGens> activeTimer;
  LifeCountdown<maxCellActiveStreakGens> streakTimer;

  LifeStableState lastTest;

  unsigned timeSincePropagate;

  unsigned currentGen;

  bool hasInteracted;
  unsigned interactionStart;
  unsigned recoveredTime;

  // Share of the search tree below this node (1 at the root, halved at each
  // branch), for the progress estimate.
  double weight = 1.0;
  unsigned depth = 0;   // number of branches above this node

  SearchParams *params;
  std::vector<Solution> *allSolutions;
  std::set<std::string> *seenRotors;
  LifeStableState *stableAtInteraction;

  SearchState(SearchParams &inparams, std::vector<Solution> &outsolutions, std::set<std::string> &outrotors, LifeStableState &stableAtInteraction);
  SearchState(const SearchState &) = default;
  SearchState &operator=(const SearchState &) = default;

  LifeState ForcedInactiveCells(
      const FrontierGeneration &gen, const LifeState &everActive,
      const LifeCountdown<maxCellActiveWindowGens> &activeTimer,
      const LifeCountdown<maxCellActiveStreakGens> &streakTimer) const;


  LifeState ForcedUnchangingCells(
      const FrontierGeneration &gen, const LifeState &everActive,
      const LifeCountdown<maxCellActiveWindowGens> &activeTimer,
      const LifeCountdown<maxCellActiveStreakGens> &streakTimer) const;


  bool UpdateActive(FrontierGeneration &generation,
                    LifeCountdown<maxCellActiveWindowGens> &activeTimer,
                    LifeCountdown<maxCellActiveStreakGens> &streakTimer);
  std::pair<bool, bool> SetForced(FrontierGeneration &generation);

  std::pair<bool, bool> TestActive(FrontierGeneration &generation);

  std::tuple<bool, bool> PopulateFrontier();
  bool UpdateFrontierStrip(unsigned column);
  bool CalculateFrontier();
  bool RefineFrontier();
  std::pair<bool, bool> TryAdvance();

  // INT port: count-based OptionsFor removed (SetForced is generic now)

  std::pair<unsigned, std::pair<int, int>> ChooseBranchCell() const;

  void SearchStep();

  void RecordOscillator();

  void RecordSolution();
  void PrintSolution(const Solution &solution);

  void SanityCheck();
};

LifeState SearchState::ForcedInactiveCells(
    const FrontierGeneration &gen,
    const LifeState &everActive,
    const LifeCountdown<maxCellActiveWindowGens> &activeTimer,
    const LifeCountdown<maxCellActiveStreakGens> &streakTimer) const {
  if (gen.gen < params->minFirstActiveGen) {
    return ~LifeState();
  }

  auto activePop = gen.active.GetPop();

  if (hasInteracted && !params->reportOscillators && gen.gen > interactionStart + params->maxActiveWindowGens)
    return ~LifeState();

  if (params->maxActiveCells != -1 && activePop > (unsigned)params->maxActiveCells)
    return ~LifeState();

  LifeState result;

  if (params->maxActiveCells != -1 && activePop == (unsigned)params->maxActiveCells)
    result |= ~gen.active; // Or maybe just return

  if (params->activeBounds.first != -1 && activePop > 0) {
    result |= ~gen.active.BufferAround(params->activeBounds);
  }

  if (params->maxEverActiveCells != -1 && everActive.GetPop() == (unsigned)params->maxEverActiveCells) {
    result |= ~everActive; // Or maybe just return
  }

  if (params->everActiveBounds.first != -1) {
    result |= ~everActive.BufferAround(params->everActiveBounds);
  }

  if (params->maxComponentActiveCells != -1 && activePop > (unsigned)params->maxComponentActiveCells) {
    for (auto &c : gen.active.Components()) {
      auto componentPop = c.GetPop();
      if(componentPop > (unsigned)params->maxComponentActiveCells)
        return ~LifeState();
      if(componentPop == (unsigned)params->maxComponentActiveCells)
        result |= ~gen.active & c.BigZOI();
    }
  }

  if (params->maxComponentEverActiveCells != -1 && everActive.GetPop() > (unsigned)params->maxComponentEverActiveCells) {
    for (auto &c : everActive.Components()) {
      auto componentPop = c.GetPop();
      if(componentPop > (unsigned)params->maxComponentEverActiveCells)
        return ~LifeState();
      if(componentPop == (unsigned)params->maxComponentEverActiveCells)
        result |= ~c & c.BigZOI();
    }
  }

  if (params->componentEverActiveBounds.first != -1) {
    for (auto &c : everActive.Components()) {
      auto wh = c.WidthHeight();
      if (wh.first > params->componentEverActiveBounds.first ||
          wh.second > params->componentEverActiveBounds.second)
        return ~LifeState();

      result |= ~c.BufferAround(params->componentEverActiveBounds) & c.BigZOI();
    }
  }

  if (params->maxCellActiveWindowGens != -1 &&
      hasInteracted &&
      gen.gen > interactionStart + (unsigned)params->maxCellActiveWindowGens)
    result |= activeTimer.finished;

  if (params->maxCellActiveStreakGens != -1 &&
      hasInteracted &&
      gen.gen > interactionStart + (unsigned)params->maxCellActiveStreakGens)
    result |= streakTimer.finished;

  if (params->maxCellStationaryDistance != -1) {
    LifeState unchanging = ~(gen.changes | (gen.next.unknown & ~gen.next.unknownStable));
    result |= unchanging.MatchLive(LifeState::NZOIAround({0, 0}, params->maxCellStationaryDistance));
  }

  return result;
}

LifeState SearchState::ForcedUnchangingCells(
    const FrontierGeneration &gen,
    const LifeState &everActive,
    const LifeCountdown<maxCellActiveWindowGens> &activeTimer,
    const LifeCountdown<maxCellActiveStreakGens> &streakTimer)
    const {
  LifeState result;
  if (params->maxChanges != -1) {
    unsigned changesPop = gen.changes.GetPop();
    if (changesPop > (unsigned)params->maxChanges)
      return ~LifeState();
    if (changesPop == (unsigned)params->maxChanges) {
      result |= ~gen.changes;
    }
  }

  if (params->maxComponentChanges != -1) {
    for (auto &c : gen.changes.Components()) {
      unsigned changesPop = c.GetPop();
      if (changesPop > (unsigned)params->maxComponentChanges)
        return ~LifeState();
      if (changesPop == (unsigned)params->maxComponentChanges) {
        result |= ~gen.changes & c.BigZOI();
      }
    }
  }

  if (params->changesBounds.first != -1) {
    result |= ~gen.changes.BufferAround(params->changesBounds);
  }

  if (params->componentChangesBounds.first != -1) {
    for (auto &c : gen.changes.Components()) {
      auto wh = c.WidthHeight();
      if (wh.first > params->componentChangesBounds.first ||
          wh.second > params->componentChangesBounds.second)
        return ~LifeState();

      result |= ~c.BufferAround(params->componentChangesBounds) & c.BigZOI();
    }
  }

  if (params->hasStator)
    result |= params->stator;

  return result;
}

bool SearchState::UpdateActive(FrontierGeneration &generation,
                               LifeCountdown<maxCellActiveWindowGens> &activeTimer,
                               LifeCountdown<maxCellActiveStreakGens> &streakTimer) {
  generation.active = generation.next.ActiveComparedTo(stable) & stable.dead0 & ~params->exempt;
  generation.changes = generation.next.ChangesComparedTo(generation.state) & stable.dead0 & ~params->exempt;

  everActive |= generation.active;

  generation.forcedInactive =
      ForcedInactiveCells(generation, everActive, activeTimer, streakTimer) &
      ~params->exempt;

  if (!(generation.active & generation.forcedInactive).IsEmpty())
    return false;

  generation.forcedUnchanging =
      ForcedUnchangingCells(generation, everActive, activeTimer, streakTimer) &
      ~params->exempt;

  if (!(generation.changes & generation.forcedUnchanging).IsEmpty())
    return false;

  return true;
}

// INT port: generic version. For a frontier cell, the current values of its
// 3x3 are known cells, or stable unknowns (cells that still equal the
// unknown stable background). Each allowed transition is checked by
// enumerating those stable unknowns: transitions with no consistent
// assignment are removed, and stable cells with the same value in every
// consistent assignment are forced. Cells with an active (non-stable)
// unknown in their 3x3 are skipped.
std::pair<bool, bool> SearchState::SetForced(FrontierGeneration &generation) {
  bool anyChanges = false;
  LifeState remainingCells = generation.frontierCells;
  for (auto cell = remainingCells.FirstOn(); cell != std::make_pair(-1, -1);
       remainingCells.Erase(cell), cell = remainingCells.FirstOn()) {

    generation.state.TransferStable(stable);
    generation.next.TransferStable(stable);

    auto allowedTransitions = generation.AllowedTransitions(stable, cell);
    if (allowedTransitions == Transition::IMPOSSIBLE)
      return {false, false};

    const LifeUnknownState &cur = generation.state;
    unsigned x = cell.first, y = cell.second;
    unsigned on9 = StableIndex9(cur.state, x, y);
    unsigned unk9 = StableIndex9(cur.unknown, x, y);
    unsigned unkStable9 = StableIndex9(cur.unknownStable, x, y);
    unsigned stableUnk9 = StableIndex9(stable.unknown, x, y);
    if ((unk9 & ~(unkStable9 & stableUnk9)) != 0)
      continue; // an active unknown nearby: nothing to deduce here

    auto pairOf = [](Transition t) -> int { // prev<<1 | next, or -1 for STABLE
      switch (t) {
      case Transition::OFF_TO_OFF: return 0;
      case Transition::OFF_TO_ON: return 1;
      case Transition::ON_TO_OFF: return 2;
      case Transition::ON_TO_ON: return 3;
      default: return -1;
      }
    };

    unsigned andMask = 511, orMask = 0;
    bool anyAssignment = false;
    Transition realised = Transition::IMPOSSIBLE;
    unsigned sub = 0;
    while (true) {
      unsigned xbits = on9 | sub;
      unsigned prev = (xbits >> 8) & 1, next = gIntRule.next[xbits];
      int pr = (int)(prev << 1 | next);
      bool ok = false;
      for (Transition t : {Transition::OFF_TO_OFF, Transition::OFF_TO_ON, Transition::ON_TO_OFF,
                           Transition::ON_TO_ON, Transition::STABLE_TO_STABLE}) {
        if ((allowedTransitions & t) != t) continue;
        int want = pairOf(t);
        bool match = want >= 0 ? want == pr : prev == next;
        if (match) { ok = true; realised |= t; }
      }
      if (ok) { anyAssignment = true; andMask &= xbits; orMask |= xbits; }
      if (sub == unk9) break;
      sub = (sub - unk9) & unk9;
    }
    if (!anyAssignment)
      return {false, false};
    allowedTransitions &= realised;

    unsigned forcedOn = andMask & unk9, forcedOff = ~orMask & unk9 & 511u;
    for (unsigned b = 0; b < 9; b++) {
      if (!((forcedOn | forcedOff) >> b & 1)) continue;
      std::pair<int, int> c = {(int)((x + N + kPosDX[b]) % N), (int)((y + 64 + kPosDY[b]) & 63)};
      if (forcedOn >> b & 1) stable.SetOn(c); else stable.SetOff(c);
      anyChanges = true;
    }
    if (forcedOn | forcedOff) {
      generation.state.TransferStable(stable);
      generation.next.TransferStable(stable);
    }

    allowedTransitions = TransitionSimplify(allowedTransitions);

    // Force the transition to occur if it is the only one allowed
    if (TransitionIsSingleton(allowedTransitions)) {
      auto transition = allowedTransitions;
      generation.SetTransition(cell, transition);
      generation.frontierCells.Erase(cell);
      if (generation.state.TransitionIsActive(cell, transition))
        everActive.Set(cell);
    }
  }

  return {true, anyChanges};
}

std::tuple<bool, bool> SearchState::PopulateFrontier() {
  bool anyChanges = false;

  frontier.state.TransferStable(stable);

  LifeUnknownState lookahead = frontier.state;
  unsigned gen = currentGen;

  auto lookaheadActiveTimer = activeTimer;
  auto lookaheadStreakTimer = streakTimer;

  for (unsigned i = 0; i < maxFrontierGens; i++) {
    gen++;

    FrontierGeneration generation;

    generation.state = lookahead;
    generation.next = lookahead.StepMaintaining(stable);
    generation.gen = gen;

    bool updateresult = UpdateActive(generation, lookaheadActiveTimer, lookaheadStreakTimer);
    if (!updateresult)
      return {false, false};

    if (params->maxCellActiveWindowGens != -1) {
      lookaheadActiveTimer.Start(generation.active);
      lookaheadActiveTimer.Tick();
    }
    if (params->maxCellActiveStreakGens != -1) {
      lookaheadStreakTimer.Reset(~generation.active);
      lookaheadStreakTimer.Start(generation.active);
      lookaheadStreakTimer.Tick();
    }

    LifeState prevUnknownActive = generation.state.unknown & ~generation.state.unknownStable;
    LifeState becomeUnknown = (generation.next.unknown & ~generation.next.unknownStable) & ~prevUnknownActive;

    generation.frontierCells = becomeUnknown & ~prevUnknownActive.ZOI();

    if (i == 0)
      frontier = generation;

    auto [result, someForced] = SetForced(generation);
    if (!result)
      return {false, false};

    anyChanges = anyChanges || someForced;

    bool isInert = ((generation.state.state ^ generation.next.state) &
                   ~generation.state.unknown & ~generation.next.unknown)
                      .IsEmpty() ||
                  (stable.dead0 & ~params->exempt & ~generation.next.unknown).IsEmpty();

    if (isInert)
      break;

    lookahead = generation.next;
  }
  return {true, anyChanges};
}

bool SearchState::UpdateFrontierStrip(unsigned column) {
  stable.SynchroniseStateKnown();

  frontier.state.TransferStable(stable);

  auto stripState = frontier.state.state.GetStrip<4>(column);
  auto stripUnknown = frontier.state.unknown.GetStrip<4>(column);
  auto stripUnknownStable = frontier.state.unknownStable.GetStrip<4>(column);

  std::tie(stripState, stripUnknown, stripUnknownStable) = frontier.state.StepMaintainingStrip(stable, column);
  frontier.next.state.SetStrip<4>(column, stripState);
  frontier.next.unknown.SetStrip<4>(column, stripUnknown);
  frontier.next.unknownStable.SetStrip<4>(column, stripUnknownStable);

  bool updateresult = UpdateActive(frontier, activeTimer, streakTimer);
  if (!updateresult)
    return false;

  return true;
}

std::pair<bool, bool> SearchState::TryAdvance() {
  bool didAdvance = false;
  bool done = false;
  while (!done) {
    if (!(frontier.next.unknown & ~frontier.next.unknownStable).IsEmpty())
      break;

    didAdvance = true;

    frontier.state = frontier.next;
    frontier.next = frontier.next.StepMaintaining(stable);
    currentGen += 1;

    LifeState active = frontier.state.ActiveComparedTo(stable) & stable.dead0 & ~params->exempt;
    everActive |= active;

    if (params->maxCellActiveWindowGens != -1) {
      activeTimer.Start(active);
      activeTimer.Tick();
    }

    if (params->maxCellActiveStreakGens != -1) {
      streakTimer.Reset(~active);
      streakTimer.Start(active);
      streakTimer.Tick();
    }

    // INT port: branching no longer marks the interaction start
    if (!hasInteracted && !active.IsEmpty()) {
      hasInteracted = true;
      interactionStart = currentGen - 1;
      *stableAtInteraction = stable;
    }

    if (hasInteracted && params->requireEscape) {
      // Exact prune: nothing anywhere differs from the stable state (and no cell is
      // undetermined), so the rest of the evolution is the still life itself (no B0):
      // nothing can escape. Checked on the whole grid, not just the catalyst zone.
      const LifeUnknownState &fs = frontier.state;
      bool undetermined = !(fs.unknown & ~fs.unknownStable).IsEmpty() ||
                          !(~fs.unknown & stable.unknown).IsEmpty();
      if (!undetermined && fs.ActiveComparedTo(stable).IsEmpty())
        return {false, false};
    }

    if (hasInteracted) {
      bool isRecovered = active.IsEmpty();

      if (isRecovered)
        recoveredTime++;
      else
        recoveredTime = 0;

      if (isRecovered && recoveredTime == params->minStableInterval) {
        if(!params->reportOscillators)
          RecordSolution();
        if(!params->continueAfterSuccess)
          return {false, false};
      }

      if (currentGen > interactionStart + params->maxActiveWindowGens) {
        if (params->reportOscillators)
          RecordOscillator();

        return {false, false};
      }
    } else {
      if (currentGen >= params->maxFirstActiveGen)
        return {false, false};
    }
  }

  return {true, didAdvance};
}

bool SearchState::CalculateFrontier() {

  unsigned rounds = 0;

  bool anyChanges = true;
  while (anyChanges) {
    anyChanges = false;

    rounds++;
    if (rounds > maxCalculateRounds)
      break;

    auto [consistent, someChanges] = PopulateFrontier();
    if (!consistent) {
      return false;
    }
    anyChanges = anyChanges || someChanges;

    // This is more important now than in old Barrister: we otherwise
    // spend a fair bit of time searching uncompletable parts of the
    // search space

    LifeState toTest = stable.Vulnerable() & stable.Differences(lastTest).ZOI();
    lastTest = stable;
    auto propagateResult = stable.TestUnknowns(toTest);

    if (!propagateResult.consistent) {
      return false;
    }
    anyChanges = anyChanges || propagateResult.changed;

    propagateResult = stable.Propagate();
    if (!propagateResult.consistent) {
      return false;
    }
    anyChanges = anyChanges || propagateResult.changed;

    stable.SanityCheck();
  }

  if (params->hasForbidden) {
    for(auto &f : params->forbiddens) {
      bool allKnown = (f.mask & stable.unknown).IsEmpty();
      if (!allKnown)
        continue;
      bool matches = ((stable.state ^ f.state) & f.mask).IsEmpty();
      if (allKnown && matches)
        return false;
    }
  }

  frontier.state.TransferStable(stable);
  frontier.next.TransferStable(stable);

  auto [consistent, didAdvance] = TryAdvance();
  if (!consistent)
    return false;

  if (didAdvance) {
    // We have to start over
    [[clang::musttail]]
    return CalculateFrontier();
  }

  return true;
}

bool SearchState::RefineFrontier() {
  stable.SynchroniseStateKnown();

  frontier.state.TransferStable(stable);
  frontier.next.TransferStable(stable);

  bool updateresult = UpdateActive(frontier, activeTimer, streakTimer);
  if (!updateresult) {
    return false;
  }

  auto [consistent, changed] = SetForced(frontier);
  if (!consistent) {
    return false;
  }

  return true;
}

std::pair<unsigned, std::pair<int, int>> SearchState::ChooseBranchCell() const {
  // unsigned stopGen = std::min(frontier.size, maxBranchingGens);

  // Prefer unknown stable cells?
  // Prefer lower numbers of possible transitions?
  // Prefer cells where all options are active?
  // for (unsigned i = 0; i < generations.size(); i++) {
  //   std::pair<int, int> branchCell = generations[i].frontierCells.FirstOn();
  //   if (branchCell.first != -1) {
  //     auto allowedTransitions = ::AllowedTransitions(generations[i], stable, branchCell);
  //     if ((allowedTransitions & Transition::STABLE_TO_STABLE) != Transition::STABLE_TO_STABLE)
  //     // if(TransitionCount(allowedTransitions) <= 3)
  //       return {i, branchCell};
  //   }
  // }
  // for (unsigned i = frontier.start; i < frontier.start + stopGen; i++) {
  //   auto &g = frontier.generations[i];
  //   auto branchCell = g.frontierCells.FirstOn();
  //   if (branchCell.first != -1)
  //     return {i, branchCell};
  // }

  LifeState remainingCells = frontier.frontierCells;
  for (auto cell = remainingCells.FirstOn(); cell != std::make_pair(-1, -1);
       remainingCells.Erase(cell), cell = remainingCells.FirstOn()) {
    auto allowedTransitions = frontier.AllowedTransitions(stable, cell);
    if(TransitionCount(allowedTransitions) <= 2)
      return {0, cell};
  }

  auto branchCell = frontier.frontierCells.FirstOn();
  if (branchCell.first != -1)
    return {0, branchCell};

  return {0, {-1, -1}};
}

struct NodeCount { unsigned long n=0, adv=0; ~NodeCount(){ if (getenv("INTDBG")) fprintf(stderr, "nodes %lu\n", n);} };
inline NodeCount gNodes;

// Progress: a node's weight counts as done when its subtree is finished, so the
// total reaches 1 when the search ends. Treats both halves of every branch as
// the same size, so it is only a rough guide (it tends to jump).
static double gProgressDone = 0;
static auto gStartTime = std::chrono::steady_clock::now();
static auto gLastProgress = std::chrono::steady_clock::now();
static size_t gSolutionsSoFar = 0;
static char gPath[33] = "";   // the current node's first 32 branch choices
void WriteResultsFile(const std::string &path, std::vector<Solution> &pats, const std::string &note);
static void ResetProgress() {
  gProgressDone = 0;
  gPath[0] = 0;
  gStartTime = gLastProgress = std::chrono::steady_clock::now();
}
static void MaybePrintProgress() {
  auto now = std::chrono::steady_clock::now();
  if (now - gLastProgress < std::chrono::seconds(10))
    return;
  gLastProgress = now;
  double secs = std::chrono::duration<double>(now - gStartTime).count();
  fprintf(stderr, "progress: %.4f%% (rough), path %s, %lu nodes, %zu solutions, %.0f s\n",
          100.0 * gProgressDone, gPath, gNodes.n, gSolutionsSoFar, secs);
}
struct ProgressGuard {   // adds the node's weight minus what its children took
  double &passed; double weight;
  ~ProgressGuard() { gProgressDone += weight - passed; }
};

void SearchState::SearchStep() {
  gNodes.n++;
  double passed = 0;
  ProgressGuard guard{passed, weight};
  if ((gNodes.n & 4095) == 0)
    MaybePrintProgress();
#ifdef DEBUG
  if (params->hasOracle) {
    if (!stable.CompatibleWith(params->oracle))
      return;
  }
#endif

  if(frontier.frontierCells.IsEmpty() || timeSincePropagate >= maxBranchFastCount){
    bool consistent = CalculateFrontier();
    if (!consistent)
      return;
    timeSincePropagate = 0;

    stable.SanityCheck();
    // SanityCheck();
  } else {
    timeSincePropagate++;

    bool consistent = RefineFrontier();
    if (!consistent)
      return;

    if (frontier.frontierCells.IsEmpty()) {
      bool consistent = CalculateFrontier();
      if (!consistent)
        return;
      timeSincePropagate = 0;
    }
  }

  auto [i, branchCell] = ChooseBranchCell();
  if (branchCell.first == -1) {
    if (getenv("INTDBG")) std::cerr << "PRUNE no-branch-cell gen " << currentGen << std::endl;
    return; // INT port: nothing to branch on (prunes; crude)
  }

  // INT port: branch on the value of an unknown stable cell in the frontier
  // cell's 3x3 instead of on the frontier cell's transition. The transition
  // branch relied on Life count reasoning to resolve the cell's next state.
  stable.SynchroniseStateKnown(branchCell);
  frontier.state.TransferStable(stable);
  LifeState nb;
  nb.Set(branchCell);
  nb = nb.ZOI() & stable.unknown;
  auto valueCell = nb.FirstOn();
  if (valueCell.first == -1) {
    LifeState z;
    z.Set(branchCell);
    z = z.ZOI();
    if ((z & frontier.state.unknown).IsEmpty()) {
      // The stored next generation is stale (its neighbourhood became known
      // after it was computed): recompute the frontier instead of pruning.
      frontier.frontierCells = LifeState();
      passed = weight;   // the repeated step below counts this node's share
      return SearchStep();
    }
    if (getenv("INTDBG")) std::cerr << "PRUNE active-unknown-neighbour gen " << currentGen << std::endl;
    return; // INT port: active unknown nearby, no stable cell to branch on (prunes; crude)
  }

  for (int v = 0; v < 2; v++) {

    SearchState newSearch = *this;
    newSearch.weight = weight / 2;
    newSearch.depth = depth + 1;
    if (depth < 32) {   // path string: 0 = first branch (cell off), 1 = second (cell on)
      gPath[depth] = '0' + v;
      gPath[depth + 1] = 0;
    }
    if (v == 0) newSearch.stable.SetOff(valueCell); else newSearch.stable.SetOn(valueCell);
    newSearch.stable.SynchroniseStateKnown(valueCell);

    auto propagateResult = newSearch.stable.PropagateStrip(valueCell.first);
    if (!propagateResult.consistent)
      continue;

    newSearch.frontier.frontierCells = LifeState(); // force a full recompute
    passed += newSearch.weight;
    newSearch.SearchStep();
  }
}

SearchState::SearchState(SearchParams &inparams,
                         std::vector<Solution> &outsolutions,
                         std::set<std::string> &outrotors,
                         LifeStableState &inStableAtInteraction)
    : currentGen{0}, hasInteracted{false}, interactionStart{0} {
  params = &inparams;
  allSolutions = &outsolutions;
  seenRotors = &outrotors;
  stableAtInteraction = &inStableAtInteraction;

  stable = inparams.stable;
  frontier.state = inparams.startingState;
  frontier.next = frontier.state.StepMaintaining(stable);

  timeSincePropagate = 0;

  everActive = LifeState();
  activeTimer = LifeCountdown<maxCellActiveWindowGens>(params->maxCellActiveWindowGens);
  streakTimer = LifeCountdown<maxCellActiveStreakGens>(params->maxCellActiveStreakGens);

  TryAdvance();
}

void SearchState::PrintSolution(const Solution &solution) {
  std::cout << "Winner:" << std::endl;
  std::cout << "x = 0, y = 0, rule = LifeBellman" << std::endl;
  LifeState state = params->startingState.state | solution.stable.state;
  LifeState marked = solution.stable.unknown | solution.stable.state;
  LifeState startingOff = (params->stable.state & ~params->startingState.state);
  state &= ~startingOff;
  marked &= ~startingOff;
  std::cout << LifeBellmanRLEFor(state, marked) << std::endl;

  switch (solution.completionResult) {
  case CompletionResult::COMPLETED:
    std::cout << "Completed:" << std::endl;
    std::cout << solution.state.RLE() << std::endl;
  break;
  case CompletionResult::INCONSISTENT:
    std::cout << "Completion Failed: Inconsistent" << std::endl;
    std::cout << LifeState().RLE() << std::endl;
  break;
  case CompletionResult::TIMEOUT:
    std::cout << "Completion Failed: Timeout" << std::endl;
    std::cout << LifeState().RLE() << std::endl;
  break;
  }
}

void SearchState::RecordOscillator() {
  unsigned period = DeterminePeriod(frontier.state, stable);
  if (period >= params->reportOscillatorsMinPeriod) {
    std::cout << "Oscillating! Period: " << period << std::endl;

    if(!(everActive.ZOI() & stable.unknown).IsEmpty()) {
      auto [result, completed] = stable.CompleteStable(
          params->stabiliseResultsTimeout, params->minimiseResults);
      if(!completed.IsEmpty()) {
        stable.SetOn(completed);
        stable.SetOff(~completed);
        frontier.state.TransferStable(stable);
      }
    }

    for(auto &r : GetSeparatedRotorDesc(frontier.state, stable, period)) {
      auto rotorDesc = r.ToString();
      if (seenRotors->contains(rotorDesc))
        std::cout << "Known Rotor: " << rotorDesc << std::endl;
      else {
        seenRotors->insert(rotorDesc);
        std::cout << "New Rotor: " << rotorDesc << std::endl;
        RecordSolution();
      }
    }
  }
}

void SearchState::RecordSolution() {
  Solution solution;
  solution.stable = stable;
  solution.interactionStable = *stableAtInteraction;
  solution.interactionGen = interactionStart;
  solution.recoveryGen = currentGen - params->minStableInterval + 1;

  if (params->stabiliseResults) {
    std::tie(solution.completionResult, solution.completed) = stable.CompleteStable(params->stabiliseResultsTimeout, params->minimiseResults);
  }

  LifeState startingActive = params->startingState.state & ~params->stable.state;
  LifeState startingStableOff = params->stable.state & ~params->startingState.state;

  solution.state = (stable.state | startingActive | solution.completed) & ~startingStableOff;

  // INT port: exact re-check, since the still-life solver and completion
  // are Life's. Only meaningful when the completion succeeded.
  if (!params->metasearch && params->stabiliseResults && solution.completionResult == CompletionResult::COMPLETED) {
    LifeState still = (stable.state | solution.completed) & ~startingStableOff;
    LifeState stepped = still;
    stepped.Step();
    bool isStill = stepped == still;

    LifeState cur = solution.state;
    unsigned T = solution.recoveryGen + params->minStableInterval;
    for (unsigned g = 0; g < T; g++)
      cur.Step();
    LifeState catalyst = still & ~(params->startingState.state & ~params->stable.state);
    bool recovered = ((cur ^ catalyst) & catalyst.ZOI()).IsEmpty();
    std::cout << "INT check: still-life=" << (isStill ? "yes" : "NO")
              << " catalyst-region-restored-at-gen-" << T << "=" << (recovered ? "yes" : "NO")
              << std::endl;
  }

  allSolutions->push_back(solution);
  gSolutionsSoFar++;

  if (!params->metasearch) {
    PrintSolution(solution);
    WriteResultsFile(params->outputFile, *allSolutions, "search still running (written after each new one)");
  }
}

void SearchState::SanityCheck() {
#ifdef DEBUG
  frontier.state.SanityCheck(stable);
  assert((stable.state & stable.unknown).IsEmpty());

  LifeStableState copy = stable;
  auto result = copy.Propagate();
  assert(result.consistent);
  assert(copy == stable); // The stable state should be fully propagated

  LifeStableState copycopy = copy;
  copycopy.Propagate();
  assert(copycopy == copy);

  LifeUnknownState startingStable = {stable.state, stable.unknown, stable.unknown};
  LifeUnknownState steppedStable = startingStable.StepMaintaining(stable);

  assert(startingStable == steppedStable);
#endif
}

void PrintSummary(std::vector<Solution> &pats, std::ostream &out) {
  out << "x = 0, y = 0, rule = " << gIntRule.name << std::endl;
  for (unsigned i = 0; i < pats.size(); i += 8) {
    std::vector<Solution> rowSolutions = std::vector<Solution>(pats.begin() + i, pats.begin() + std::min((unsigned)pats.size(), i + 8));
    std::vector<LifeState> row;
    for (auto &s : rowSolutions) {
      row.push_back(s.state);
    }
    out << RowRLE(row) << std::endl;
  }
}

void PrintSummary(std::vector<Solution> &pats) {
  PrintSummary(pats, std::cout);
}

// Rewrites the results file (via a temporary file, so it is never half-written).
void WriteResultsFile(const std::string &path, std::vector<Solution> &pats, const std::string &note) {
  if (path.empty())
    return;
  std::string tmp = path + ".tmp";
  {
    std::ofstream f(tmp);
    f << "#C " << pats.size() << " solution(s), " << note << std::endl;
    PrintSummary(pats, f);
  }
  std::rename(tmp.c_str(), path.c_str());
}

bool PassesFilters(const SearchParams &params, const Solution &solution) {
  LifeUnknownState state = params.startingState;
  state.TransferStable(solution.stable);

  unsigned filterTime = 0;
  for (auto &f : params.filters) {
    filterTime = std::max(filterTime, f.gen);
  }

  std::vector<bool> filterPassed(params.filters.size(), false);
  for (unsigned i = 0; i < filterTime; i++) {
    state = state.StepMaintaining(solution.stable);

    if (i < solution.interactionGen)
      continue;

    for (unsigned fi = 0; fi < params.filters.size(); fi++) {
      auto &f = params.filters[fi];

      if (!(state.unknown & f.mask).IsEmpty())
        break;

      bool shouldCheck = f.type == FilterType::EVER || (f.type == FilterType::EXACT && i+1 == f.gen);
      if (shouldCheck && ((state.state ^ f.state) & f.mask).IsEmpty())
        filterPassed[fi] = true;
    }
  }

  for (bool b : filterPassed) {
    if (!b)
      return false;
  }
  return true;
}

std::vector<Solution> TrimSolutions(SearchParams &params, std::vector<Solution> &solutions) {
  unsigned maxGen = params.maxFirstActiveGen + params.maxActiveWindowGens + params.minStableInterval;

  // A list of hashes per generation
  std::vector<std::vector<std::pair<uint64_t, Solution>>> hashes(maxGen);
  for (auto &s : solutions) {
    LifeStableState clearedStable = s.stable.ClearUnmodified();
    LifeUnknownState state = params.startingState;
    state.TransferStable(clearedStable);

    LifeState stator = clearedStable.dead0 & ~clearedStable.unknown;

    // Fast forward to when the catalyst has just recovered
    for (unsigned i = 0; i < s.recoveryGen; i++) {
      state = state.StepMaintaining(clearedStable);
      stator &= ~state.ActiveComparedTo(clearedStable);
    }

    s.stator = stator;

    // Add hashes to the list until the catalyst is destroyed/interacted with a second time
    for (unsigned i = s.recoveryGen; i < maxGen; i++) {
      bool isRecovered = ((clearedStable.state ^ state.state) & clearedStable.dead0 & ~params.exempt).IsEmpty();
      if (!isRecovered) {
        break;
      }

      LifeState perturbed = state.state & ~clearedStable.dead0;

      uint64_t hash = perturbed.GetHash();
      bool hashSeen = false;
      for (auto &[oldhash, oldsolution] : hashes[i]) {
        if (hash == oldhash) {
          if (s < oldsolution) {
            oldsolution = s;
          }
          hashSeen = true;
        }
      }
      if (!hashSeen) {
        hashes[i].push_back({hash, s});
      }

      state = state.StepMaintaining(clearedStable);
    }
  }

  std::map<Solution, unsigned> counts;
  for (auto &g : hashes) {
    for (auto &[hash, s] : g) {
      counts[s]++;
    }
  }

  std::vector<Solution> results;
  for (auto &[s, count] : counts) {
    if (count >= params.minTrimHashes)
      results.push_back(s);
  }

  return results;
}

void MetaSearchStep(unsigned round, std::vector<Solution> &allSolutions, SearchParams &params) {
  std::vector<Solution> roundSolutions;
  std::set<std::string> seenRotors;
  LifeStableState stableAtInteraction;

  std::cerr << "Depth: " << round << std::endl;
  std::cerr << "x = 0, y = 0, rule = LifeBellman" << std::endl;
  std::cerr << LifeBellmanRLEFor(params.stable.state | params.startingState.state, params.stable.unknown | params.stable.state) << std::endl;

  SearchState search(params, roundSolutions, seenRotors, stableAtInteraction);
  ResetProgress();
  search.SearchStep();

  auto trimmed = TrimSolutions(params, roundSolutions);

  std::vector<Solution> filtered;
  for (auto &s : trimmed) {
    if (PassesFilters(params, s))
      filtered.push_back(s);
  }
  allSolutions.insert(allSolutions.end(), filtered.begin(), filtered.end());

  if (filtered.size() > 0) {
    std::cout << "Winner!" << std::endl;
    PrintSummary(filtered);
    WriteResultsFile(params.outputFile, allSolutions, "metasearch still running (written after each round)");
  }

  if (round == params.metasearchRounds) {
    return;
  }

  for (auto &s : trimmed) {
    SearchParams newParams = params;
    newParams.minFirstActiveGen = s.interactionGen;
    newParams.stable = s.interactionStable.Graft(s.stable);
    newParams.startingState.TransferStable(newParams.stable);
    newParams.exempt |= newParams.stable.dead0 & ~s.stator;
    newParams.stator |= s.stator;
    if (round == 1) {
      newParams.minFirstActiveGen = std::max(s.interactionGen, params.minMetaFirstActiveGen);
      newParams.maxFirstActiveGen = params.maxMetaFirstActiveGen;
    }

    MetaSearchStep(round + 1, allSolutions, newParams);
  }
}

void MetaSearch(SearchParams &params) {
  std::vector<Solution> allSolutions;
  MetaSearchStep(1, allSolutions, params);
  std::cout << "Summary!" << std::endl;
  PrintSummary(allSolutions);
  WriteResultsFile(params.outputFile, allSolutions, "metasearch complete");
}

int main(int, char *argv[]) {
  auto toml = toml::parse(argv[1]);
  SearchParams params = SearchParams::FromToml(toml);
  if (params.outputFile.empty()) {   // default: next to the input, <input>-results.rle
    std::string base = argv[1];
    if (base.size() > 5 && base.compare(base.size() - 5, 5, ".toml") == 0)
      base.resize(base.size() - 5);
    params.outputFile = base + "-results.rle";
  }
  std::cerr << "Results are written to " << params.outputFile << " as they are found" << std::endl;

  if (params.maxCellActiveWindowGens != -1 &&
  (unsigned)params.maxCellActiveWindowGens > maxCellActiveWindowGens) {
    std::cout << "max-cell-active-window is higher than allowed by the hardcoded value!" << std::endl; exit(1);
  }
  if (params.maxCellActiveStreakGens != -1 &&
  (unsigned)params.maxCellActiveStreakGens > maxCellActiveStreakGens) {
    std::cout << "max-cell-active-streak is higher than allowed by the hardcoded value!" << std::endl; exit(1);
  }

  if (params.metasearch) {
    MetaSearch(params);
  } else {
    std::vector<Solution> allSolutions;
    LifeStableState stableAtInteraction;

    std::set<std::string> seenRotors;
    if(params.reportOscillators) {
      for (auto r : ReadRotors(params.knownrotorsFile))
      seenRotors.insert(r);
    }

    SearchState search(params, allSolutions, seenRotors, stableAtInteraction);
    ResetProgress();
    search.SearchStep();
    fprintf(stderr, "progress: %.4f%%, %lu nodes, %zu solutions, search complete\n", 100.0 * gProgressDone, gNodes.n, allSolutions.size());
    WriteResultsFile(params.outputFile, allSolutions, "search complete");

    if (params.printSummary) {
      std::cout << "All solutions:" << std::endl;
      PrintSummary(allSolutions);

      if(params.trimResults) {
        std::cout << "Unique perturbations:" << std::endl;
        auto trimmed = TrimSolutions(params, allSolutions);
        PrintSummary(trimmed);
        allSolutions = trimmed;
      }

      if(!params.filters.empty()) {
        std::vector<Solution> filtered;
        for (auto &s : allSolutions) {
          if (PassesFilters(params, s))
            filtered.push_back(s);
        }
        std::cout << "Filtered:" << std::endl;
        PrintSummary(filtered);
      }
    }
  }
}
