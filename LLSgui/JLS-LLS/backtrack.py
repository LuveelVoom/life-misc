#!/usr/bin/env python3
"""Runs a JLS status file in native/llsbt, the lifesrc-style multithreaded backtracker.

Usage: backtrack.py <file.jdf> [-n SOLUTIONS] [--threads N] [-o FILE] [LLS arguments...]

The problem comes from LLS's JLS import and constraints (through builder.py, which is
LLS-GUI's builder.py). LLS then hands the rule to the backtracker as a table rather than as
clauses. Needs a fully specified rule without B0.
"""

import argparse
import os
import shlex
import subprocess
import sys
import time

import builder

ROOT = os.path.dirname(os.path.realpath(__file__))


def fail(message):
    print(message)
    sys.exit(1)


parser = argparse.ArgumentParser()
parser.add_argument("jdf")
parser.add_argument("-n", "--number_of_solutions", type=int, default=1)
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("-o", "--output_file_name")
args, lls_args = parser.parse_known_args()

with open(args.jdf) as f:
    pattern, jdf_args = builder.import_jdf(f.read())
lls_args = jdf_args + lls_args
if "--layer_max_rotor" in lls_args:
    fail("The backtracker can't do JLS active cells per layer limits: use a SAT solver")

start = time.time()
try:
    problem = builder.build(pattern, shlex.join(lls_args), record_transitions=True)
except builder.UnsatInPreprocessing:
    print("Unsatisfiable (in preprocessing)")
    sys.exit(0)
except builder.BuildError as e:
    fail(str(e))
if problem.search_pattern.rule["B0c"] != "0":
    fail("The backtracker hasn't been tested with B0 rules: use a SAT solver")

print("Variables: %d, clauses: %d, preprocessing %.1f s" % (
    problem.number_of_variables, problem.number_of_clauses, time.time() - start))
print("Solving with %d thread%s...\n" % (args.threads, "" if args.threads == 1 else "s"), flush=True)

base = os.path.realpath(args.jdf)
paths = [base + suffix for suffix in (".cnf", ".watch", ".trans")]
try:
    try:
        problem.write_backtrack_files(*paths)
    except builder.BuildError as e:
        fail(str(e))
    process = subprocess.Popen(
        [os.path.join(ROOT, "native", "llsbt")] + paths
        + ["--threads", str(args.threads), "--solutions", str(args.number_of_solutions), "--interval", "3600000"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)  # EOF on stdin would mean stop
    found = 0
    for line in process.stdout:
        kind, _, rest = line.rstrip("\n").partition(" ")
        if kind == "SOL":
            _, elapsed, values = rest.split(" ", 2)
            rle = problem.decode(values)["rle"]
            found += 1
            print(rle + "\n", flush=True)
            if args.output_file_name:
                with open(args.output_file_name, "a") as f:
                    f.write(rle + "\n")
        elif kind == "UNSAT":
            print("Search exhausted" if found else "Unsatisfiable", flush=True)
        elif kind == "STOPPED":
            print("Stopped", flush=True)
        elif kind == "DONE":
            print("Solutions: %s, solver time: %s s" % (rest.split()[1], rest.split()[0]), flush=True)
            break
    process.kill()
finally:
    for path in paths:
        if os.path.exists(path):
            os.remove(path)
