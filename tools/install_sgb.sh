#!/bin/sh
# Download and build Stanford GraphBase locally into ext/sgb/.
# This avoids the need for distribution-packaged libgb-dev.
# Works on Linux and macOS.

set -e

mkdir -p ext
cd ext

if [ -f sgb/libgb.a ]; then
    echo "ext/sgb/libgb.a already built; skipping download."
    exit 0
fi

if [ ! -f sgb.tar.gz ]; then
    echo "Downloading Stanford GraphBase..."
    curl -fsSL -o sgb.tar.gz https://cs.stanford.edu/~knuth/sgb.tar.gz
fi

if [ ! -d sgb ]; then
    mkdir sgb
    tar xf sgb.tar.gz -C sgb --strip-components=0
fi

cd sgb
echo "Building libgb.a..."
make lib > /dev/null
ls -la libgb.a
echo ""
echo "Done. The library lives at ext/sgb/libgb.a; headers at ext/sgb/."
