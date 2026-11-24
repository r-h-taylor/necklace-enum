#!/bin/sh
# Richard H. Taylor, necklace-enum, doi:10.5281/zenodo.24330095
# A short tour of necklace-enum from the command line.
# pip install necklace-enum ase
# sh necklace-try.sh
#
# The first run after installation takes a few extra seconds while Numba compiles the
# generator to machine code. Later runs do not require the load time. 
set -e
if [ -x /usr/bin/time ]; then
  if [ "$(uname)" = "Darwin" ]; then T="/usr/bin/time -l"; else T="/usr/bin/time -v"; fi
else
  T=""   # no /usr/bin/time: peak memory is not shown
fi

echo "== 1. Reproduce Hart & Forcade (2008), Table V: binary fcc, up to label exchange, n = 2-23"
necklace-enum count fcc 2 2-23 --label-exchange | tail -4

echo
echo "== 2. Beyond that table, by complete generation: n = 24-28 (Table 3 of the paper)"
necklace-enum count fcc 2 24-28 --label-exchange

echo
echo "== 3. A large fixed composition: fcc 24:8 at n = 32 (22.7 million structures)."
echo "   enumlib needs ~50 GB for this problem. To compare, watch the peak memory reported below."
$T necklace-enum count fcc Cu,Au 32 --composition 24,8 2>&1 | awk '
  /n= /                       { print }
  /peak memory footprint/     { printf "peak memory: %.0f MB\n", $1 / 1e6 }        # macOS: bytes
  /Maximum resident set size/ { printf "peak memory: %.0f MB\n", $NF / 1e3 }       # Linux: kbytes
'

echo
echo "== 4. List every one of them, as integer arrays in Python (no structure is stored)"
python3 - <<'PY'
import time, necklace_enum as ne
P = ne.fcc(["Cu", "Au"])
t = time.perf_counter()
n = sum(labs.shape[0] for _, labs in ne.enumerate_labelings(P, 32, composition=(24, 8)))
print(f"listed {n:,d} structures in {time.perf_counter() - t:.1f} s")
PY

echo
echo "== 5. Write structures for other codes: hcp Mg5Cd3 in 4-cell supercells, as extended XYZ"
necklace-enum list hcp Mg,Cd 4 --composition 5,3 --format extxyz -o MgCd_n4.xyz
grep -c Lattice MgCd_n4.xyz | xargs echo "structures written to MgCd_n4.xyz:"

echo
echo "== 6. Drop-in for enumlib: necklace-enum reads struct_enum.in and writes struct_enum.out"
cat > struct_enum.in <<'IN'
fcc
bulk
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
 2 -nary case
 1
 0.0 0.0 0.0 0/1
 1 10
0.10000000E-06
full list of labelings
IN
necklace-enum enumlib struct_enum.in
echo "(enumlib's makeStr.py reads this struct_enum.out unchanged)"
