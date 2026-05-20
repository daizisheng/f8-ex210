# Makefile for the f8-ex210 toolkit. Supports Linux and macOS.
#
# Stanford GraphBase (libgb) is needed by make_knight and the two
# patched DYNAHAM binaries. The Makefile prefers a locally-built copy
# at ext/sgb/ (produced by ./tools/install_sgb.sh) and falls back to
# the system libgb-dev package if it is installed.

CC      ?= cc
SRCDIR   = src

# ---- OS-specific flags ------------------------------------------------------

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Large code model is needed on Linux x86_64 because of DYNAHAM's
# large static arrays.
ifeq ($(UNAME_S),Linux)
  ifeq ($(UNAME_M),x86_64)
    LARGE = -mcmodel=large
  endif
endif

# OpenMP: gcc/clang on Linux ship with it. Apple's clang doesn't, so
# users with `brew install libomp` should set, e.g.,
#   make OMP="-Xpreprocessor -fopenmp -lomp"
OMP ?= -fopenmp

# Pick up the locally-built libgb if present; otherwise rely on the
# system installation.
ifneq ($(wildcard ext/sgb/libgb.a),)
  SGB_CFLAGS  = -I ext/sgb
  SGB_LDFLAGS = ext/sgb/libgb.a
else
  SGB_CFLAGS  =
  SGB_LDFLAGS = -lgb
endif

CFLAGS_OPT = -O3 $(OMP)
CFLAGS_DYN = -O2 $(LARGE) $(SGB_CFLAGS)
LDFLAGS_DYN = $(SGB_LDFLAGS)

# ---- Targets ---------------------------------------------------------------

.PHONY: all clean

all: simulate_count divtest make_knight dynaham_modp_dumpT dynahamp_modp_dumpT verify_bm

verify_bm: $(SRCDIR)/verify_bm.c
	$(CC) -O3 $< -o $@

simulate_count: $(SRCDIR)/simulate_count.c
	$(CC) $(CFLAGS_OPT) $< -o $@

divtest: $(SRCDIR)/divtest.c
	$(CC) -O3 $< -o $@

make_knight: $(SRCDIR)/make_knight.c
	$(CC) -O2 $(SGB_CFLAGS) $< $(LDFLAGS_DYN) -o $@

dynaham_modp_dumpT: $(SRCDIR)/dynaham_modp_dumpT.c
	$(CC) $(CFLAGS_DYN) $< $(LDFLAGS_DYN) -o $@

dynahamp_modp_dumpT: $(SRCDIR)/dynahamp_modp_dumpT.c
	$(CC) $(CFLAGS_DYN) $< $(LDFLAGS_DYN) -o $@

clean:
	rm -f simulate_count divtest make_knight verify_bm
	rm -f dynaham_modp_dumpT dynahamp_modp_dumpT
	rm -f *.gb *.log Q5_tilde.txt Q5_closed.txt Q5_open.txt
