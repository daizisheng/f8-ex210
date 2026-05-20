#!/bin/sh
# Download Knuth's dynaham.w and dynahamp.w from his website.
# Saved under ext/.

set -e

mkdir -p ext
cd ext

for f in dynaham.w dynahamp.w; do
    if [ -f "$f" ]; then
        echo "already have $f"
    else
        echo "fetching $f ..."
        curl -fsS -o "$f" "https://cs.stanford.edu/~knuth/programs/$f"
    fi
done

ls -la dynaham.w dynahamp.w
