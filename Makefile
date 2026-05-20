CC      = cc
CFLAGS  = -O3 -fopenmp
LDFLAGS =

BINDIR  = .
SRCDIR  = src

# Stanford GraphBase library (libgb) and FLINT
LDFLAGS_DYNAHAM = -lgb
LDFLAGS_DIVTEST = -lflint

# Larger memory model needed for DYNAHAM's static arrays.
CFLAGS_DYNAHAM  = -O2 -mcmodel=large

.PHONY: all clean

all: simulate_count divtest dynaham_modp_dumpT dynahamp_modp_dumpT make_knight

simulate_count: $(SRCDIR)/simulate_count.c
	$(CC) $(CFLAGS) $< -o $@

divtest: $(SRCDIR)/divtest.c
	$(CC) -O3 $< $(LDFLAGS_DIVTEST) -o $@

dynaham_modp_dumpT: $(SRCDIR)/dynaham_modp_dumpT.c
	$(CC) $(CFLAGS_DYNAHAM) $< $(LDFLAGS_DYNAHAM) -o $@

dynahamp_modp_dumpT: $(SRCDIR)/dynahamp_modp_dumpT.c
	$(CC) $(CFLAGS_DYNAHAM) $< $(LDFLAGS_DYNAHAM) -o $@

make_knight: $(SRCDIR)/make_knight.c
	$(CC) -O2 $< $(LDFLAGS_DYNAHAM) -o $@

clean:
	rm -f simulate_count divtest dynaham_modp_dumpT dynahamp_modp_dumpT make_knight
	rm -f *.gb *.log
