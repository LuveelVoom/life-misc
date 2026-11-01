"""Builds LLS search problems from GUI requests, using LLS itself as a library.

The argument handling mirrors the `lls` command-line script, so any LLS
argument (symmetries, population limits, partial rules, backgrounds, ...)
can be passed through from the GUI.
"""

import argparse
import contextlib
import itertools
import io
import os
import re
import shlex
import sys
import threading

APP_DIR = os.path.dirname(os.path.realpath(__file__))
LLS_DIR = os.path.join(APP_DIR, "lls")
sys.path.insert(0, LLS_DIR)

import settings  # noqa: E402
import src.formatting  # noqa: E402
import src.literal_manipulation  # noqa: E402
import src.logging  # noqa: E402
import src.rules  # noqa: E402
import src.search_patterns  # noqa: E402
from src.literal_manipulation import negate, neighbours_from_coordinates, variable_from_literal  # noqa: E402
from src.SearchPattern import SearchPattern  # noqa: E402
from src.UnsatInPreprocessing import UnsatInPreprocessing  # noqa: E402

# LLS keeps logging state in module globals, so only build one problem at a time
_build_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Faster transition encoding (LLS method 2, used for every rule except plain Life)
#
# LLS builds all 1024 clauses per cell through generic helpers and then drops
# those satisfied by constant cells or constant rule transitions. This produces
# the same clause set, but only enumerates neighbourhoods consistent with the
# constants and writes the DIMACS strings directly.
# ---------------------------------------------------------------------------

_transition_of = {}


def _transition_table():
    if not _transition_of:
        for alive in itertools.product([True, False], repeat=8):
            _transition_of[alive] = src.rules.transition_from_cells(alive)
    return _transition_of


def _force_transition_method_2(self, grid, x, y, t):
    table = _transition_table()
    clauses = self.clauses
    dimacs_of = clauses.dimacs_literal_from_variable
    duration = len(grid)
    cell = grid[t][y][x]
    predecessor = grid[(t - 1) % duration][y][x]
    # LLS's own method 2 takes neighbours from self.grid even for the background's cells,
    # unlike methods 0 and 1; use the grid being evolved, as they do
    neighbours = neighbours_from_coordinates(grid, x, y, t, background_grid=self.background_grid)

    def info(literal):
        """(variable number, positive DIMACS literal string, negative one) for a non-constant literal."""
        negated = literal[0] == "-"
        variable = literal[1:] if negated else literal
        number = dimacs_of.get(variable)
        if number is None:
            clauses.number_of_variables += 1
            number = str(clauses.number_of_variables)
            dimacs_of[variable] = number
        positive, negative = ("-" + number, number) if negated else (number, "-" + number)
        return number, positive, negative

    def options(literal):
        # (alive?, clause literal or None). "alive" puts the negated literal in the clause.
        if literal == "0":
            return [(False, None)]
        if literal == "1":
            return [(True, None)]
        number, positive, negative = info(literal)
        return [(True, (number, negative)), (False, (number, positive))]

    def side(literal):
        # clause literal for "literal" (None if constant false; "SAT" if constant true)
        if literal == "0":
            return None
        if literal == "1":
            return "SAT"
        number, positive, _ = info(literal)
        return (number, positive)

    def negation(literal):
        return "1" if literal == "0" else "0" if literal == "1" else (literal[1:] if literal[0] == "-" else "-" + literal)

    cell_true, cell_false = side(cell), side(negation(cell))
    rule_cache = {}
    clause_set = clauses.clause_set

    def emit(body, extra):
        seen = {}
        for item in body + extra:
            if item is None:
                continue
            if item == "SAT":
                return
            number, text = item
            previous = seen.get(number)
            if previous is None:
                seen[number] = text
            elif previous != text:
                return  # tautology
        literals = sorted(seen.values())
        literals.append("0\n")
        clause_set.add(" ".join(literals))

    position_options = [options(predecessor)] + [options(n) for n in neighbours]
    for combination in itertools.product(*position_options):
        alive = tuple(o[0] for o in combination[1:])
        key = ("S" if combination[0][0] else "B") + table[alive]
        transition = self.rule[key]
        if key not in rule_cache:
            rule_cache[key] = (side(negation(transition)), side(transition))
        not_transition, yes_transition = rule_cache[key]
        body = [o[1] for o in combination]
        # transition holds -> cell is alive; transition fails -> cell is dead
        emit(body, [not_transition, cell_true])
        emit(body, [yes_transition, cell_false])


# ---------------------------------------------------------------------------
# Method 3 (from the newer LLS): one clause per prime implicate of the rule's
# transition function on (cell, predecessor, 8 neighbours), so unit propagation
# draws every conclusion a single cell's neighbourhood allows. Same clauses as
# LLS's method 3, but constants are simplified while writing.
# ---------------------------------------------------------------------------

_primes_of = {}


def _make_prime(n, fn):
    """Prime implicates of fn over n bits, as (values, don't-care mask) pairs.
    Copied from the newer LLS's src/primer.py (which needs Python 3.12)."""
    t = []
    for i in range(1 << n):
        if not fn([(i >> j & 1) == 1 for j in range(n)]):
            t.append([i, 0, 0])
    for c in t:
        for d in t:
            diff = c[0] ^ d[0]
            if diff and diff & (diff - 1) == 0:
                c[2] |= diff
    for v in range(n):
        s = []
        j = 0
        for c in t:
            if not c[2]:
                yield c[0], c[1]
            else:
                if c[2] >> (v + 1):
                    s.append(c)
                if not c[0] >> v & 1 and c[2] >> v & 1:
                    while t[j][0] != c[0] | 1 << v or t[j][1] != c[1]:
                        j += 1
                    s.append([c[0], c[1] | 1 << v, c[2] & t[j][2] & ~(1 << v)])
        t = s
    for c in t:
        if not c[2]:
            yield c[0], c[1]


def _rule_primes(rule):
    key = tuple(sorted(rule.items()))
    if key not in _primes_of:
        def holds(bits):
            transition = ("S" if bits[1] else "B") + src.rules.transition_from_cells(bits[2:])
            return rule[transition] == ("1" if bits[0] else "0")

        # (a, b): the clause over the positions not in mask b that excludes their values in a
        _primes_of[key] = list(_make_prime(10, holds))
    return _primes_of[key]


def _force_transition_method_3(self, grid, x, y, t):
    clauses = self.clauses
    dimacs_of = clauses.dimacs_literal_from_variable
    clause_set = clauses.clause_set
    duration = len(grid)
    inputs = [grid[t][y][x], grid[(t - 1) % duration][y][x]]
    inputs += neighbours_from_coordinates(grid, x, y, t, background_grid=self.background_grid)
    # per input: None if constant, else (variable number, DIMACS literal true when it is alive)
    resolved = []
    for literal in inputs:
        if literal in ("0", "1"):
            resolved.append(literal == "1")
            continue
        negated = literal[0] == "-"
        variable = literal[1:] if negated else literal
        number = dimacs_of.get(variable)
        if number is None:
            clauses.number_of_variables += 1
            number = str(clauses.number_of_variables)
            dimacs_of[variable] = number
        resolved.append((number, ("-" + number) if negated else number, number if negated else "-" + number))
    for a, b in _rule_primes(self.rule):
        seen = {}
        for i in range(10):
            if b >> i & 1:
                continue
            alive = a >> i & 1
            item = resolved[i]
            if item is True or item is False:
                if item != bool(alive):
                    break  # satisfied by a constant
                continue
            number, when_alive, when_dead = item
            text = when_dead if alive else when_alive
            previous = seen.get(number)
            if previous is None:
                seen[number] = text
            elif previous != text:
                break  # tautology
        else:
            literals = sorted(seen.values())
            literals.append("0\n")
            clause_set.add(" ".join(literals))


_original_force_transition = SearchPattern.force_transition

# When not None, transitions are recorded here (for the backtracking engine) instead of
# being encoded as clauses. Only touched while holding _build_lock.
_recorded_transitions = None


def _register(clauses, literal):
    """Give a literal's variable a DIMACS number even if no clause mentions it."""
    if literal in ("0", "1"):
        return
    variable, _ = variable_from_literal(literal)
    if variable not in clauses.dimacs_literal_from_variable:
        clauses.number_of_variables += 1
        clauses.dimacs_literal_from_variable[variable] = str(clauses.number_of_variables)


def _force_transition(self, grid, x, y, t, method):
    if _recorded_transitions is not None:
        duration = len(grid)
        cell = grid[t][y][x]
        predecessor = grid[(t - 1) % duration][y][x]
        neighbours = neighbours_from_coordinates(grid, x, y, t, background_grid=self.background_grid)
        literals = [predecessor] + list(neighbours) + [cell]
        for literal in literals:
            _register(self.clauses, literal)
        _recorded_transitions.append(literals)
        return
    if method == 3:
        if all(value in "01" for value in self.rule.values()):
            return _force_transition_method_3(self, grid, x, y, t)
        method = 2  # partial rules: prime implicates need a fixed rule
    if method == 2:
        return _force_transition_method_2(self, grid, x, y, t)
    return _original_force_transition(self, grid, x, y, t, method)


SearchPattern.force_transition = _force_transition

# Default transition encoding when -M isn't given (LLS's own default: 1 for Life, 2 otherwise)
DEFAULT_METHOD = None


def _force_evolution(self, method=None):
    """LLS's force_evolution, also accepting method 3."""
    if method is None:
        method = DEFAULT_METHOD
    if method is None:
        method = settings.life_encoding_method if src.rules.rulestring_from_rule(self.rule) == "B3/S23" else 2
    if method not in range(4):
        raise BuildError("Method not found: " + str(method))
    if method < 2 and src.rules.rulestring_from_rule(self.rule) != "B3/S23":
        raise BuildError("Rules other than Life can only use method 2 or 3")
    src.logging.log("Enforcing evolution rule...", 1)
    src.logging.log("Method: " + str(method))
    starting_number_of_clauses = len(self.clauses.clause_set)
    for t, generation in enumerate(self.grid):
        if t > 0:
            for y, row in enumerate(generation):
                for x, cell in enumerate(row):
                    if not self.ignore_transition[t][y][x]:
                        self.force_transition(self.grid, x, y, t, method)
    for t, generation in enumerate(self.background_grid):
        for y, row in enumerate(generation):
            for x, cell in enumerate(row):
                if not self.background_ignore_transition[t][y][x]:
                    self.force_transition(self.background_grid, x, y, t, method)
    src.logging.log("Number of clauses used: " + str(len(self.clauses.clause_set) - starting_number_of_clauses))
    src.logging.log("Done\n", -1)


SearchPattern.force_evolution = _force_evolution


class BuildError(Exception):
    pass


class Aborted(Exception):
    """A build was abandoned because a newer request superseded it."""


def make_parser():
    parser = argparse.ArgumentParser(prog="lls", add_help=False, exit_on_error=False)
    parser.add_argument('-s', '--symmetry', nargs="*", action="append", default=[])
    parser.add_argument('-a', '--asymmetry', nargs="*", action="append", default=[])
    parser.add_argument('-c', '--force_change', nargs="*", action="append", type=int, default=[])
    parser.add_argument("-p", "--population", "--pop", action='append', nargs="*", default=[])
    parser.add_argument('-M', '--method', type=int, default=None)
    parser.add_argument("-r", '--rule', type=str, default=None)
    parser.add_argument('--background', "--bg", nargs=1, default=None)
    parser.add_argument('--background_offset', "--bgos", nargs=3, type=int, default=None)
    parser.add_argument('--max_change', type=int, default=None)
    parser.add_argument('--max_decay', type=int, default=None)
    parser.add_argument('--max_growth', type=int, default=None)
    parser.add_argument('--layer_max_rotor', type=int, default=None)
    parser.add_argument('--layer_shape', default=None)
    parser.add_argument('--layer_center', type=int, nargs=2, default=None)
    parser.add_argument("--layer_rotor_include_set_cells", action="store_true")
    parser.add_argument('--polyomino', type=int, nargs="?", default=None, const=-1)
    return parser


def list_backgrounds():
    return sorted(f for f in os.listdir(os.path.join(LLS_DIR, "backgrounds")) if not f.startswith("."))


def split_args(text):
    """Split an LLS-style argument string, keeping backslashes (used by RE\\ symmetries)."""
    lexer = shlex.shlex(text, posix=True)
    lexer.whitespace_split = True
    lexer.escape = ""
    return list(lexer)


def parse_symmetries(argument_list):
    """Same logic as the lls script: turn -s/-a arguments into [transformation, x, y, period]."""
    result = []
    max_period = 0
    for arguments in argument_list:
        transformations = []
        period = x_translate = y_translate = None
        for argument in arguments:
            re_match = re.match("(^.*?)([0-9]*$)", argument)
            letters = re_match.group(1).upper()
            number = re_match.group(2)
            if letters == "P":
                if period is not None:
                    raise BuildError("Can only have one period per symmetry")
                period = int(number)
                max_period = max(max_period, period)
            elif letters == "X":
                x_translate = int(number) if number else 0
            elif letters == "Y":
                y_translate = int(number) if number else 0
            elif letters == "X-":
                x_translate = -int(number)
            elif letters == "Y-":
                y_translate = -int(number)
            else:
                if transformations:
                    raise BuildError("Can only impose one transformation per symmetry: " + " ".join(arguments))
                transformation = letters + str(number)
                named = {
                    "C1": ["RO0"], "C2": ["RO2"], "C4": ["RO1"],
                    "D2-": ["RE-"], "D2\\": ["RE\\"], "D2|": ["RE|"], "D2/": ["RE/"],
                    "D4+": ["RE-", "RE|"], "D4X": ["RE/", "RE\\"], "D8": ["RE-", "RE\\"],
                }
                if letters == "RO":
                    transformations.append("RO" + str(int(number or 0) % 4))
                elif transformation in named:
                    transformations += named[transformation]
                elif transformation in ["RE-", "RE/", "RE|", "RE\\"]:
                    transformations.append(transformation)
                else:
                    raise BuildError('Symmetry argument "' + transformation + '" not recognized')
        if not transformations:
            transformations = ["RO0"]
        for transformation in transformations:
            result.append([transformation, x_translate or 0, y_translate or 0, period or 0])
    return result, max_period


def parse_population(arguments_list):
    at_most, at_least, exactly = [], [], []
    for arguments in arguments_list:
        if len(arguments) == 0:
            at_least.append([[0], 1])
            continue
        constraint = arguments[0]
        times = [int(t) for t in arguments[1:]] if len(arguments) > 1 else [0]
        re_match = re.match("(^.*?)([0-9]+$)", constraint)
        if not re_match:
            raise BuildError('Population constraint "' + constraint + '" not recognized')
        sign = re_match.group(1)
        amount = int(re_match.group(2))
        if sign in ["", "="]:
            exactly.append([times, amount])
        elif sign in ["<=", "=<"]:
            at_most.append([times, amount])
        elif sign in [">=", "=>"]:
            at_least.append([times, amount])
        elif sign == "<":
            at_most.append([times, amount - 1])
        elif sign == ">":
            at_least.append([times, amount + 1])
        else:
            raise BuildError('Population constraint "' + constraint + '" not recognized')
    return at_most, at_least, exactly


class Problem:
    """A preprocessed LLS search, ready for a SAT solver."""

    def __init__(self, search_pattern, determined, show_background, log_text, transitions=None):
        self.search_pattern = search_pattern
        self.transitions = transitions
        self.determined = determined
        self.show_background = show_background
        self.log_text = log_text
        clauses = search_pattern.clauses
        dimacs_of = clauses.dimacs_literal_from_variable

        # Every variable that appears in the grid, background or rule gets watched,
        # so the GUI can display it live and solutions can be decoded from it.
        watched = []
        index_of = {}

        def code_of(literal):
            if literal == "0":
                return 0
            if literal == "1":
                return 1
            variable, negated = variable_from_literal(literal)
            if variable not in dimacs_of:
                return -1  # unconstrained cell
            dimacs = int(dimacs_of[variable])
            if dimacs not in index_of:
                index_of[dimacs] = len(watched)
                watched.append(dimacs)
            return 2 + 2 * index_of[dimacs] + (1 if negated else 0)

        grid = search_pattern.grid
        self.cell_codes = [[[code_of(c) for c in row] for row in gen] for gen in grid]
        self.background_codes = [[[code_of(c) for c in row] for row in gen]
                                 for gen in search_pattern.background_grid]
        self.rule_codes = {t: code_of(l) for t, l in sorted(search_pattern.rule.items())}
        self.watched = watched

        # Literals that distinguish solutions. LLS's force_distinct uses only generation 0
        # when the pattern is determined by it; using every generation gives the same
        # solutions then, and lets the backtracker decide cells in any order.
        block = []
        for gen in grid:
            for row in gen:
                block += [c for c in row]
        for gen in search_pattern.background_grid:
            for row in gen:
                block += [c for c in row]
        block += list(search_pattern.rule.values())
        block_dimacs = set()
        for literal in block:
            if literal in ["0", "1"]:
                continue
            variable, _ = variable_from_literal(literal)
            if variable in dimacs_of:
                block_dimacs.add(int(dimacs_of[variable]))
        self.block = sorted(block_dimacs)

        self.width = len(grid[0][0]) - 2
        self.height = len(grid[0]) - 2
        self.duration = len(grid)
        self.number_of_variables = clauses.number_of_variables
        self.number_of_clauses = len(clauses.clause_set)
        self.number_of_cells = search_pattern.number_of_cells()

    def dimacs(self):
        return self.search_pattern.clauses.make_string()

    def write_files(self, cnf_path, watch_path):
        with open(cnf_path, "w") as f:
            f.write(self.dimacs())
        with open(watch_path, "w") as f:
            f.write(" ".join(map(str, self.watched)) + "\n")
            f.write(" ".join(map(str, self.block)) + "\n")

    def write_backtrack_files(self, cnf_path, watch_path, transitions_path, order="auto", origin=(0, 0)):
        """Input for native/llsbt: the clauses, plus each cell transition and the rule table."""
        rule = self.search_pattern.rule
        bits = []
        for letter in "BS":
            for mask in range(256):
                alive = tuple(bool(mask >> i & 1) for i in range(8))
                value = rule[letter + src.rules.transition_from_cells(alive)]
                if value not in ("0", "1"):
                    raise BuildError("The backtracking engine needs a fully specified rule; "
                                     "use the SAT engine for partial rules.")
                bits.append(value)
        dimacs_of = self.search_pattern.clauses.dimacs_literal_from_variable

        def dimacs(literal):
            if literal == "0":
                return "F"
            if literal == "1":
                return "T"
            variable, negated = variable_from_literal(literal)
            return ("-" if negated else "") + dimacs_of[variable]

        self.write_files(cnf_path, watch_path)
        with open(transitions_path, "w") as f:
            f.write("".join(bits) + "\n")
            f.write(" ".join(map(str, self.decision_order(order, origin))) + "\n")
            for literals in self.transitions:
                f.write(" ".join(dimacs(l) for l in literals) + "\n")

    def decision_order(self, order, origin=(0, 0)):
        """Cell literals meaning "dead", in the order the backtracker should decide them.

        Orders: "auto", "rows", "columns" decide one generation at a time. "jls-<shape>"
        follows JLS: all generations of a cell together, cells ordered by distance from
        `origin` (inner coordinates), measured as columns |dx|, rows |dy|, box, diamond or
        circle. Variables that distinguish solutions come first (the backtracker relies
        on this to skip duplicate solutions)."""
        grid = self.search_pattern.grid
        dimacs_of = self.search_pattern.clauses.dimacs_literal_from_variable
        duration, height, width = len(grid), len(grid[0]), len(grid[0][0])
        inner = [(x, y) for y in range(1, height - 1) for x in range(1, width - 1)]
        if order == "auto":
            # scan along the shorter side of the field first (lifesrc does the same)
            order = "rows" if width <= height else "columns"
        if order.startswith("jls-"):
            ox, oy = origin[0] + 1, origin[1] + 1
            distance = {
                "columns": lambda dx, dy: abs(dx), "rows": lambda dx, dy: abs(dy),
                "box": lambda dx, dy: max(abs(dx), abs(dy)), "diamond": lambda dx, dy: abs(dx) + abs(dy),
                "circle": lambda dx, dy: dx * dx + dy * dy,
            }.get(order[4:])
            if distance is None:
                raise BuildError("Unknown search order " + order)
            inner.sort(key=lambda p: (distance(p[0] - ox, p[1] - oy), (p[0] - ox) ** 2 + (p[1] - oy) ** 2,
                                      p[0], p[1]))
            coordinates = [(t, y, x) for x, y in inner for t in range(duration)]
        elif order == "columns":
            coordinates = [(t, y, x) for t in range(duration) for x in range(1, width - 1)
                           for y in range(1, height - 1)]
        else:
            coordinates = [(t, y, x) for t in range(duration) for y in range(1, height - 1)
                           for x in range(1, width - 1)]
        block = set(self.block)
        first, rest = [], []
        for t, y, x in coordinates:
            literal = grid[t][y][x]
            if literal in ("0", "1"):
                continue
            variable, negated = variable_from_literal(literal)
            number = int(dimacs_of[variable])
            dead = number if negated else -number
            (first if number in block else rest).append(dead)
        # remaining distinguishing variables (border/background/rule) before anything else
        listed = {abs(l) for l in first}
        first += [-v for v in self.block if v not in listed]
        return first + rest

    def summary(self):
        return {
            "width": self.width,
            "height": self.height,
            "duration": self.duration,
            "variables": self.number_of_variables,
            "clauses": self.number_of_clauses,
            "cells": self.number_of_cells,
            "watched": len(self.watched),
            "determined": self.determined,
            "cell_codes": self.cell_codes,
            "background_codes": self.background_codes,
            "rule_codes": self.rule_codes,
            "rulestring": src.rules.rulestring_from_rule(self.search_pattern.rule),
            "grid_text": self.preview_text(),
        }

    def preview_text(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.search_pattern.make_string(pattern_output_format="csv", show_background=True)

    def decode(self, values):
        """Turn a string of watched-variable values into an LLS solution."""
        solution = set()
        for dimacs, ch in zip(self.watched, values):
            solution.add(str(dimacs) if ch in "1t" else "-" + str(dimacs))
        with _build_lock, contextlib.redirect_stdout(io.StringIO()):
            src.logging.verbosity_level = 0
            sol = self.search_pattern.substitute_solution(solution)
            rle = sol.make_string(pattern_output_format="rle", determined=self.determined,
                                  show_background=self.show_background)
            csv = sol.make_string(pattern_output_format="csv", determined=self.determined,
                                  show_background=self.show_background)
            rulestring = src.rules.rulestring_from_rule(sol.rule)
        cells = [[[int(c) for c in row[1:-1]] for row in gen[1:-1]] for gen in sol.grid]
        return {"rle": rle, "csv": csv, "rule": rulestring, "cells": cells}


def build(pattern_text, args_text="", background_text=None, verbosity=2, should_abort=None,
          record_transitions=False):
    """Preprocess a search. Returns a Problem, or raises BuildError / UnsatInPreprocessing.

    With record_transitions, the rule is not encoded as clauses; the Problem lists each
    cell transition instead (for the backtracking engine)."""
    global _recorded_transitions
    parser = make_parser()
    try:
        args = parser.parse_args(split_args(args_text))
    except (argparse.ArgumentError, SystemExit, ValueError) as e:
        raise BuildError("Could not parse LLS arguments: " + str(e))

    # Wait for the lock, but give up if the caller no longer needs this build
    while not _build_lock.acquire(timeout=0.05):
        if should_abort and should_abort():
            raise Aborted()
    try:
        if should_abort and should_abort():
            raise Aborted()
        out = io.StringIO()
        src.logging.verbosity_level = verbosity
        src.logging.indent_level = 0
        _recorded_transitions = [] if record_transitions else None
        try:
            with contextlib.redirect_stdout(out):
                problem = _build(args, pattern_text, background_text)
            problem.transitions = _recorded_transitions
        except AssertionError as e:
            raise BuildError("LLS rejected the search: " + (str(e) or "assertion failed") +
                             "\n\n" + out.getvalue()[-3000:])
        except UnsatInPreprocessing:
            raise
        except BuildError:
            raise
        except Exception as e:
            raise BuildError("LLS error: " + type(e).__name__ + ": " + str(e) + "\n\n" + out.getvalue()[-3000:])
        problem.log_text = out.getvalue()
        return problem
    finally:
        _recorded_transitions = None
        _build_lock.release()


def _build(args, pattern_text, background_text):
    force_change = []
    for generations in args.force_change:
        if len(generations) == 0:
            force_change.append([0, 1])
        elif len(generations) == 2:
            force_change.append(generations)
        else:
            raise BuildError("Wrong number of arguments to -c")

    population_at_most, population_at_least, population_exactly = parse_population(args.population)
    symmetries, _ = parse_symmetries(args.symmetry)
    asymmetries, _ = parse_symmetries(args.asymmetry)

    rule = src.rules.rule_from_rulestring(args.rule or "B3/S23")
    if args.rule is not None and args.rule.strip()[0] == "{":
        for transition, literal in rule.items():
            if literal not in ["0", "1"]:
                variable, negated = variable_from_literal(literal)
                rule[transition] = negate("user_input_" + variable, negated)

    if background_text is None:
        name = args.background[0] if args.background else "vacuum"
        path = os.path.join(LLS_DIR, "backgrounds", os.path.basename(name))
        if not os.path.isfile(path):
            raise BuildError('Background "' + name + '" not found')
        with open(path) as f:
            background_text = f.read()
    background_grid, background_ignore = src.search_patterns.search_pattern_from_string(background_text)
    if args.background_offset:
        src.literal_manipulation.offset_background(background_grid, *args.background_offset)
        src.literal_manipulation.offset_background(background_ignore, *args.background_offset)

    grid, ignore_transition = src.search_patterns.search_pattern_from_string(pattern_text)

    search_pattern = SearchPattern(
        grid,
        ignore_transition=ignore_transition,
        background_grid=background_grid,
        background_ignore_transition=background_ignore,
        rule=rule,
    )

    search_pattern.standardise_variables_names()
    for symmetry in symmetries:
        search_pattern.force_symmetry(symmetry)
    search_pattern.remove_redundancies()
    search_pattern.standardise_variables_names()

    for asymmetry in asymmetries:
        search_pattern.force_asymmetry(asymmetry)
    for constraint in population_at_most:
        search_pattern.force_population_at_most(constraint)
    for constraint in population_at_least:
        search_pattern.force_population_at_least(constraint)
    for constraint in population_exactly:
        search_pattern.force_population_exactly(constraint)
    if args.max_change is not None:
        search_pattern.force_max_change(args.max_change)
    if args.max_decay is not None:
        search_pattern.force_max_decay(args.max_decay)
    if args.max_growth is not None:
        search_pattern.force_max_growth(args.max_growth)
    if args.layer_max_rotor is not None:
        no_subperiods = [[0] * len(grid[0][0]) for _ in grid[0]]
        search_pattern.force_layer_max_rotor(args.layer_max_rotor, args.layer_shape or "column",
                                             args.layer_center or [0, 0],
                                             args.layer_rotor_include_set_cells, no_subperiods)
    if args.polyomino is not None:
        search_pattern.force_polyomino(args.polyomino)
    for times in force_change:
        search_pattern.force_change(times)

    search_pattern.force_evolution(method=args.method)

    determined = search_pattern.deterministic()
    show_background = search_pattern.background_nontrivial()
    return Problem(search_pattern, determined, show_background, "")


def import_jdf(text):
    """Convert a JLS .jdf status file into LLS pattern text plus LLS arguments."""
    with _build_lock, contextlib.redirect_stdout(io.StringIO()):
        src.logging.verbosity_level = 0
        pattern = src.formatting.parse_jdf(text)
        lls_args, _, _ = src.formatting.get_jdf_args(text)
    return pattern, lls_args
