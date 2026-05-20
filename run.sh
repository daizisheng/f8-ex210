#!/bin/sh
# run.sh -- end-to-end pipeline for the m=5 investigation of Exercise 210.
#
# Total wall time on a 16-core Xeon Platinum 8375C: about 3 minutes.
# Set CORES below if you want to use fewer threads.

set -e

CORES=${OMP_NUM_THREADS:-16}
PRIME_FLAG=""        # set to -DMODP=1000000009ULL for the second-prime check.

echo "=== 1. Build (if not already) ==="
make all

echo ""
echo "=== 2. Generate a 5x30 knight graph ==="
./make_knight 5 30 k5x30.gb

echo ""
echo "=== 3. DYNAHAM dump runs (closed + open) ==="
./dynaham_modp_dumpT  k5x30.gb 2> Tdump_closed.log  > /dev/null
./dynahamp_modp_dumpT k5x30.gb 2> Tdump_open.log    > /dev/null
wc -l Tdump_closed.log Tdump_open.log

echo ""
echo "=== 4. Period iterator: closed sequence ==="
OMP_NUM_THREADS=$CORES OMP_PROC_BIND=close ./simulate_count \
    Tdump_closed.log 12000 40 30
mv /tmp/simulate_count.poly Q5_tilde.txt   # u-space polynomial, degree 4106

echo ""
echo "=== 5. Period iterator: open sequence ==="
OMP_NUM_THREADS=$CORES OMP_PROC_BIND=close ./simulate_count \
    Tdump_open.log 40000 40 30
mv /tmp/simulate_count.poly Q5_open.txt    # z-space polynomial, degree 39630

echo ""
echo "=== 6. Convert u-space Q to z-space (insert zeros, double indices) ==="
awk '/^[0-9]/ {print 2*$1, $2}' Q5_tilde.txt > Q5_closed.txt

echo ""
echo "=== 7. Divisibility test (Q^k divides Q+) for k = 1, 2, 3 ==="
echo "--- k = 1 ---"
./divtest Q5_closed.txt Q5_open.txt 1
echo "--- k = 2 ---"
./divtest Q5_closed.txt Q5_open.txt 2
echo "--- k = 3 ---"
./divtest Q5_closed.txt Q5_open.txt 3

echo ""
echo "Pipeline finished."
