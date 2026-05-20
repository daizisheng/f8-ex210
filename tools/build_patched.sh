#!/bin/sh
# Apply our patches to Knuth's dynaham.c / dynahamp.c (the ctangle outputs),
# producing the mod-p + T/C/W-dumping versions used by the fast path.

set -e

./tools/fetch_knuth.sh

cd ext

ctangle dynaham.w  > /dev/null
ctangle dynahamp.w > /dev/null

echo "Applying patches..."
patch -p0 -i ../patches/dynaham.patch  dynaham.c
patch -p0 -i ../patches/dynahamp.patch dynahamp.c

# rename to make their role clear
mv dynaham.c  dynaham_modp_dumpT.c
mv dynahamp.c dynahamp_modp_dumpT.c

echo ""
echo "Building..."
cc -O2 -mcmodel=large dynaham_modp_dumpT.c  -lgb -o dynaham_modp_dumpT
cc -O2 -mcmodel=large dynahamp_modp_dumpT.c -lgb -o dynahamp_modp_dumpT

cd ..

# Move resulting binaries into the top-level for run.sh.
cp ext/dynaham_modp_dumpT  .
cp ext/dynahamp_modp_dumpT .
ls -la dynaham_modp_dumpT dynahamp_modp_dumpT
