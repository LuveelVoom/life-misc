// Parses a Hensel INT string like "B2ei3aeikqr/S2-a3". Fills birth[256] and
// survive[256] (indexed by 8-bit neighbourhood mask). Returns false + err on
// failure.
static bool ParseHenselRule(const std::string &rule, bool *birth,
                            bool *survive, std::string &err) {
  std::string r;
  for (char ch : rule)
    if (!isspace((unsigned char)ch))
      r += ch;
  if (r.empty()) {
    err = "empty rule";
    return false;
  }

  bool hasBS = false;
  for (char ch : r)
    if (ch == 'B' || ch == 'b' || ch == 'S' || ch == 's')
      hasBS = true;

  // 0 = birth, 1 = survival
  int mode = -1;
  size_t i = 0;
  if (!hasBS) {
    // numeric S/B notation, e.g. 23/3
    size_t slash = r.find('/');
    if (slash == std::string::npos) {
      err = "expected something like B3/S23";
      return false;
    }
    r = "S" + r.substr(0, slash) + "/B" + r.substr(slash + 1);
  }

  bool seenB = false, seenS = false;
  while (i < r.size()) {
    char ch = r[i];
    if (ch == '/') {
      i++;
      continue;
    }
    if (ch == 'B' || ch == 'b') {
      mode = 0;
      seenB = true;
      i++;
      continue;
    }
    if (ch == 'S' || ch == 's') {
      mode = 1;
      seenS = true;
      i++;
      continue;
    }
    if (mode < 0 || ch < '0' || ch > '8') {
      err = std::string("unexpected character '") + ch + "'";
      return false;
    }
    i++;
    bool negate = false;
    if (i < r.size() && r[i] == '-') {
      negate = true;
      i++;
    }
    std::string letters;
    while (i < r.size() && isalpha((unsigned char)r[i]) && r[i] != 'B' &&
           r[i] != 'b' && r[i] != 'S' && r[i] != 's') {
      letters += (char)tolower((unsigned char)r[i]);
      i++;
    }
    bool *dest = (mode == 0) ? birth : survive;

    // validate letters
    for (char l : letters) {
      bool ok = false;
      for (int c = 0; c < NUM_HENSEL_CLASSES; c++)
        if (HENSEL_CLASSES[c].n == ch && HENSEL_CLASSES[c].l == l)
          ok = true;
      if (!ok) {
        err = std::string("letter '") + l + "' is not valid for neighbour count " + ch;
        return false;
      }
    }
    if (negate && letters.empty()) {
      err = std::string("'-' after ") + ch + " must be followed by letters";
      return false;
    }

    for (int c = 0; c < NUM_HENSEL_CLASSES; c++) {
      if (HENSEL_CLASSES[c].n != ch)
        continue;
      bool listed = letters.find(HENSEL_CLASSES[c].l) != std::string::npos;
      bool on;
      if (letters.empty())
        on = true;               // bare digit: every configuration
      else if (negate)
        on = !listed;
      else
        on = listed;
      if (on)
        for (int m = 0; m < HENSEL_CLASSES[c].cnt; m++)
          dest[HENSEL_CLASSES[c].masks[m]] = true;
    }
  }
  if (!seenB && !seenS) {
    err = "rule must contain B and/or S";
    return false;
  }
  return true;
}

static int Base64Val(char c) {
  if (c >= 'A' && c <= 'Z') return c - 'A';
  if (c >= 'a' && c <= 'z') return c - 'a' + 26;
  if (c >= '0' && c <= '9') return c - '0' + 52;
  if (c == '+') return 62;
  if (c == '/') return 63;
  return -1;
}

// Parses a Golly MAP rule (512 bits, base64). Golly's bit index is
// 256*NW + 128*N + 64*NE + 32*W + 16*C + 8*E + 4*SW + 2*S + SE, bits are read
// MSB first from each decoded byte. Fills birth/survive indexed by 8-bit mask.
static bool ParseMapRule(const std::string &rule, bool *birth, bool *survive,
                         std::string &err) {
  std::string r = rule.substr(3);
  while (!r.empty() && (r.back() == '=' || isspace((unsigned char)r.back())))
    r.pop_back();
  std::vector<int> bits;
  for (char ch : r) {
    int v = Base64Val(ch);
    if (v < 0) {
      err = "invalid base64 character in MAP rule";
      return false;
    }
    for (int k = 5; k >= 0; k--)
      bits.push_back((v >> k) & 1);
  }
  if (bits.size() < 512) {
    err = "MAP rule is too short (need 512 bits / 86 base64 characters)";
    return false;
  }
  for (int k = 0; k < 512; k++) {
    if (!bits[k])
      continue;
    int mask = ((k & 0x1e0) >> 1) | (k & 0x0f);
    if (k & 0x10)
      survive[mask] = true;
    else
      birth[mask] = true;
  }
  return true;
}
