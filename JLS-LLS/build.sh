#!/bin/sh
# Builds JavaLifeSearch_LLS.jar and native/llsbt. Needs a JDK (9 or later), a C++ compiler, and
# the JLS 1.7 download, JLS-1.7.zip (included; JLS_ZIP overrides it).
set -e
cd "$(dirname "$0")"
ZIP="${JLS_ZIP:-JLS-1.7.zip}"
[ -f "$ZIP" ] || { echo "JLS download not found at $ZIP: set JLS_ZIP"; exit 1; }
rm -rf build
mkdir -p build/jls-orig build/classes
unzip -q -o "$ZIP" 'JavaLifeSearch_1.7.jar' -d build
unzip -q -o build/JavaLifeSearch_1.7.jar 'jls/*' 'images/*' -d build/jls-orig
find build/jls-orig -name '*.class' -delete
tr -d '\r' < build/jls-orig/jls/ContentPane.java > build/ContentPane.java
mv build/ContentPane.java build/jls-orig/jls/ContentPane.java
patch -s -d build/jls-orig -p1 < ContentPane.patch
cp src/jls/LlsRunDialog.java build/jls-orig/jls/
find build/jls-orig -name '*.java' > build/sources.txt
javac -nowarn --release 9 -encoding ISO-8859-1 -d build/classes @build/sources.txt
cp -R build/jls-orig/images build/classes/
printf 'Main-Class: jls.JavaLifeSearch\n' > build/manifest.txt
jar cfm JavaLifeSearch_LLS.jar build/manifest.txt -C build/classes jls -C build/classes images
rm -rf build
c++ -O3 -std=c++17 native/llsbt.cpp -o native/llsbt -lpthread
