# transcensus-int: Apple clang + Homebrew libomp (Apple clang has no OpenMP of its own)
CC      ?= clang
LIBOMP  ?= /opt/homebrew/opt/libomp
CFLAGS  ?= -std=c99 -O3 -mcpu=native -Wall -Wextra -Wno-unused-function
OMPFLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
OMPLIBS  = -L$(LIBOMP)/lib -lomp

all: transcensus-int

transcensus-int: transcensus-int.c lib/*.c
	$(CC) $(CFLAGS) $(OMPFLAGS) transcensus-int.c -o $@ $(OMPLIBS) -lm

# Fast build fixed to one rule (about 3x faster than the generic build): make fixed RULE=B2-a3/S23
# gives transcensus-int-B2-a3_S23, which refuses any other rule. The rule's BDD is written out as C by --emit-c
NAME = $(subst /,_,$(subst +,-,$(RULE)))
fixed: transcensus-int
	@test -n "$(RULE)" || { echo "usage: make fixed RULE=<rule>"; exit 1; }
	mkdir -p build
	./transcensus-int -r '$(RULE)' --emit-c > 'build/rule-$(NAME).h'
	$(CC) $(CFLAGS) $(OMPFLAGS) -I. -DRULE_GENERATED='"build/rule-$(NAME).h"' transcensus-int.c -o 'transcensus-int-$(NAME)' $(OMPLIBS) -lm

# Single-threaded build without OpenMP
transcensus-int-serial: transcensus-int.c lib/*.c
	$(CC) $(CFLAGS) transcensus-int.c -o $@ -lm

# The unmodified original (B3/S23 only), used as the oracle in tests/
tc-orig: orig/transcensus.c orig/lib/*.c
	$(CC) -std=c99 -O3 -Wno-unused-function orig/transcensus.c -o $@ -lm

clean:
	rm -f transcensus-int transcensus-int-serial tc-orig transcensus-int-*
	rm -rf build

.PHONY: all clean fixed
