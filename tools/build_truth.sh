#!/bin/sh
# Build Knuth's dynaham (unpatched), run on small boards, print S_{5,n} ground-truth values.
# These integer values can be reduced mod p and compared with our simulate_count output.

set -e

./tools/fetch_knuth.sh

cd ext
echo "=== ctangle Knuth's bignum dynaham ==="
ctangle dynaham.w > /dev/null
echo ""

echo "=== build Knuth's bignum dynaham (closed) ==="
cc -O2 -mcmodel=large dynaham.c -lgb -o dynaham_truth
echo ""

echo "=== ctangle Knuth's bignum dynahamp ==="
ctangle dynahamp.w > /dev/null

echo "=== build Knuth's bignum dynahamp (open) ==="
cc -O2 -mcmodel=large dynahamp.c -lgb -o dynahamp_truth
echo ""

cd ..

# Generate a small 5x10 knight graph (50 vertices) -- runs in seconds
echo "=== generate 5x10 knight graph ==="
./make_knight 5 10 k5x10.gb 2>/dev/null
echo ""

echo "=== Knuth's CLOSED dynaham on 5x10 (bignum) ==="
./ext/dynaham_truth k5x10.gb 2>/dev/null | grep "Hamiltonian" | head
echo ""

echo "=== Knuth's OPEN dynahamp on 5x10 (bignum) ==="
./ext/dynahamp_truth k5x10.gb 2>/dev/null | grep "Hamiltonian" | head
echo ""

echo "These integer values are the ground truth for S_{5,n}."
echo "Compare them mod p with the first iterations of simulate_count."
