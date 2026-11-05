#!/bin/sh
# Double-click (macOS) or run from a terminal to start JLS with the LLS menu item.
cd "$(dirname "$0")"
[ -f JavaLifeSearch_LLS.jar ] || ./build.sh
exec java -jar JavaLifeSearch_LLS.jar "$@"
