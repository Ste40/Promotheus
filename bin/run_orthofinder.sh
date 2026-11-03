#!/usr/bin/env bash
set -euo pipefail
# The frozen executable prepends its own bin/ to PATH. Use a private runtime
# beside our FastME adapter, leaving the pinned image's executables untouched.
runtime="$(mktemp -d /tmp/promotheus-orthofinder.XXXXXX)"
trap 'rm -rf "$runtime"' EXIT
mkdir "$runtime/bin"
cp /OrthoFinder/orthofinder /OrthoFinder/config.json "$runtime/"
cp /adapters/fastme_deterministic.py "$runtime/bin/fastme"
chmod +x "$runtime/bin/fastme"
ln -s /OrthoFinder/bin/diamond "$runtime/bin/diamond"
ln -s /OrthoFinder/bin/mcl "$runtime/bin/mcl"
export PATH="$runtime/bin:$PATH"
"$runtime/orthofinder" "$@"
