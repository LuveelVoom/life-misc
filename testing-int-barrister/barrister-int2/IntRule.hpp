#pragma once
// INT / MAP rule support for the crude INT port of Barrister.
//
// Neighbourhood convention (matches LifeAPI here: x = column index, y = bit
// index, y grows downward): 8-bit mask NW=0x80 N=0x40 NE=0x20 W=0x10 E=0x08
// SW=0x04 S=0x02 SE=0x01; 9-bit index = centre<<8 | mask.
//
// Stepping and the still-life solver both use the exact rule tables below, so
// any INT or MAP rule without B0 is supported.

#include <cctype>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace intrule_detail {
#include "hensel_classes.inc"
#include "hensel_parse.inc"
}

struct IntRule {
  std::string name = "B3/S23";
  bool isLife = true;
  uint8_t next[512];          // centre<<8 | mask8 -> next state
  // Three-valued step: index (unk9 << 9) | on9 over 9-bit (centre<<8|mask)
  // cell sets; 0 = off, 1 = on, 2 = depends on the unknown cells.
  std::vector<uint8_t> tri;
  // Neighbour-count shadows (bit n = some neighbourhood with n live
  // neighbours does this transition), used where Barrister reasons by counts.
  unsigned offoff, offon, onoff, onon;
  // Still-life constraint at one cell: index (unk9 << 9) | on9 as above.
  // bit 31: some completion of the unknowns is stable at this cell (rule
  // applied to the 3x3 gives back the centre); bits 0-8: unknown positions
  // that are ON in every such completion; bits 9-17: unknown positions OFF in
  // every such completion.
  std::vector<uint32_t> still;
};

inline IntRule gIntRule;

inline bool IntRuleSetup(const std::string &ruleIn, std::string &err, bool unsafeSurvival = false) {
  std::string rule;
  for (char ch : ruleIn)
    if (!isspace((unsigned char)ch))
      rule += ch;

  bool birth[256] = {false}, survive[256] = {false};
  bool ok;
  if (rule.size() > 3 && (rule.compare(0, 3, "MAP") == 0 || rule.compare(0, 3, "map") == 0))
    ok = intrule_detail::ParseMapRule(rule, birth, survive, err);
  else
    ok = intrule_detail::ParseHenselRule(rule, birth, survive, err);
  if (!ok)
    return false;
  if (birth[0]) {
    err = "rules with B0 are not supported";
    return false;
  }

  IntRule &R = gIntRule;
  R.name = rule;
  R.isLife = true;
  R.offoff = R.offon = R.onoff = R.onon = 0;
  for (int m = 0; m < 256; m++) {
    int n = __builtin_popcount(m);
    R.next[m] = birth[m];
    R.next[256 | m] = survive[m];
    if (birth[m] != (n == 3) || survive[m] != (n == 2 || n == 3))
      R.isLife = false;
    (birth[m] ? R.offon : R.offoff) |= 1u << n;
    (survive[m] ? R.onon : R.onoff) |= 1u << n;
  }

  if (unsafeSurvival)
    std::cerr << "note: unsafe-allow-unrepresented-survival is no longer needed "
                 "(the still-life solver is generic) and is ignored" << std::endl;

  R.tri.assign(1u << 18, 0);
  for (unsigned unk = 0; unk < 512; unk++) {
    for (unsigned on = 0; on < 512; on++) {
      if (on & unk) continue;
      bool seen0 = false, seen1 = false;
      // enumerate subsets of unk
      unsigned sub = 0;
      while (true) {
        unsigned idx = on | sub;
        if (R.next[idx]) seen1 = true; else seen0 = true;
        if (sub == unk) break;
        sub = (sub - unk) & unk;
      }
      R.tri[(unk << 9) | on] = seen0 && seen1 ? 2 : (seen1 ? 1 : 0);
    }
  }

  R.still.assign(1u << 18, 0);
  for (unsigned unk = 0; unk < 512; unk++) {
    for (unsigned on = 0; on < 512; on++) {
      if (on & unk) continue;
      unsigned andMask = 511, orMask = 0;
      bool any = false;
      unsigned sub = 0;
      while (true) {
        unsigned x = on | sub;
        if (R.next[x] == ((x >> 8) & 1)) { any = true; andMask &= x; orMask |= x; }
        if (sub == unk) break;
        sub = (sub - unk) & unk;
      }
      uint32_t v = 0;
      if (any) v = (1u << 31) | (andMask & unk) | ((~orMask & unk & 511u) << 9);
      R.still[(unk << 9) | on] = v;
    }
  }
  return true;
}

// Build the 9-bit index for the cell at bit y of column C (L/R are the
// neighbouring columns). Torus wrap in y, as in LifeAPI.
inline unsigned IntIndex9(uint64_t L, uint64_t C, uint64_t R, unsigned y) {
  unsigned u = (y + 63) & 63, d = (y + 1) & 63;
  auto b = [](uint64_t col, unsigned yy) -> unsigned { return (col >> yy) & 1; };
  return (b(C, y) << 8) | (b(L, u) << 7) | (b(C, u) << 6) | (b(R, u) << 5) |
         (b(L, y) << 4) | (b(R, y) << 3) | (b(L, d) << 2) | (b(C, d) << 1) |
         b(R, d);
}

// Exact next state of the cells in `todo` (bits of column C).
inline uint64_t IntStepBits(uint64_t L, uint64_t C, uint64_t R, uint64_t todo) {
  uint64_t out = 0;
  while (todo) {
    unsigned y = __builtin_ctzll(todo);
    todo &= todo - 1;
    out |= (uint64_t)gIntRule.next[IntIndex9(L, C, R, y)] << y;
  }
  return out;
}
