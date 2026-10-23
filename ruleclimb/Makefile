# make        -> multi-threaded build (OpenMP). On macOS the newest Homebrew g++-N is used, since Apple clang has no OpenMP.
# make nomp   -> single-threaded build with any C++17 compiler (e.g. Apple clang)
SRC   = rulesearch.cpp
FLAGS = -O2 -march=native -std=c++17
ifeq ($(origin CXX),default)
  CXX := $(lastword g++ $(shell ls /opt/homebrew/bin/g++-[0-9]* /usr/local/bin/g++-[0-9]* 2>/dev/null | sort -t- -k2 -n))
endif

rulesearch: $(SRC)
	$(CXX) $(FLAGS) -fopenmp -o $@ $(SRC)

nomp: $(SRC)
	$(CXX) $(FLAGS) -o rulesearch $(SRC)

.PHONY: nomp
