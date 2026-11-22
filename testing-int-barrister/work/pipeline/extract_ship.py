# usage: extract_ship.py history.txt out.rle  -> exit 0 if a ship was found (after the last w9 command line)
import sys, re
L=open(sys.argv[1], errors="replace").read().split("\n")
start=max([i for i,l in enumerate(L) if "./qfind -r B3/S23-q4i -v c/7 -w 9 -s odd" in l] or [0])
ships=[]; i=start
while i < len(L):
    if L[i].startswith("Longest partial result"): break   # qfind's end-of-search partial, not a ship
    if re.match(r"^x = \d+, y = \d+, rule = ", L[i]):
        blk=[L[i]]; i+=1
        while i < len(L) and not blk[-1].rstrip().endswith("!"):
            blk.append(L[i]); i+=1
        if blk[-1].rstrip().endswith("!"): ships.append("\n".join(blk))
        continue
    i+=1
if ships:
    open(sys.argv[2],"w").write("#C qfind -r B3/S23-q4i -v c/7 -w 9 -s odd\n" + "\n\n".join(ships) + "\n")
    sys.exit(0)
sys.exit(1)
