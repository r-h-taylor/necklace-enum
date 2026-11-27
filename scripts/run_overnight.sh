#!/bin/sh
# Overnight reruns for the paper, on an otherwise idle, externally cooled machine (about 2.5 h).
#
#   cd ~/Projects/necklace-enum
#   nohup caffeinate -i sh scripts/run_overnight.sh < /dev/null > scripts/overnight_M2.out 2>&1 &
#
# Progress: tail scripts/overnight_M2.out
#
#   1. build and version record: necklace-enum, Numba, Python, enumlib commit and build flags
#   2. enumlib labdegn check, with the per-size Totdegn totals
#   3. Table 1 timings of necklace-enum with this version (arrays with copy=True, objects without
#      degeneracies); the 0.3.2 log is kept in superseded/
#   4. binary fcc at n = 32 counted up to label exchange (the "230 s" in the text)
#   5. enumlib on fcc 42:6 at n = 48, with its full output, to see how it fails
cd "$(dirname "$0")/.." || exit 1
ENUMLIB="${ENUMLIB:-$HOME/Projects/enumlib}"
ENUMX="$ENUMLIB/src/enum.x"
S=scripts
echo "started $(date)"

echo "== 1. versions and builds"
{
  date
  python3 -c "import platform, numba, numpy, necklace_enum as te; print(platform.platform()); print('Python', platform.python_version(), '| NumPy', numpy.__version__, '| Numba', numba.__version__); print('necklace-enum', te.__version__, 'from', te.__file__)"
  echo "necklace-enum repository: $(git log -1 --format='%h %cd' --date=short)"
  echo "enumlib repository: $(git -C "$ENUMLIB" log -1 --format='%h %cd %s' --date=short)"
  git -C "$ENUMLIB" status --short
  echo
  bash $S/check_enumlib_build.sh "$ENUMLIB"
} > $S/build_info_M2.txt 2>&1

echo "== 2. enumlib labdegn check"
python3 $S/check_enumlib_labdegn.py "$ENUMX" > $S/enumlib_labdegn_check_M2.txt 2>&1

echo "== 3. necklace-enum timings for Table 1"
[ -f $S/necklace_timings_M2.txt ] && mv $S/necklace_timings_M2.txt $S/superseded/necklace_timings_M2_v032.txt
python3 $S/time_necklace.py --cooldown 60 > $S/necklace_timings_M2.txt 2>&1

echo "== 4. binary fcc, n = 32, up to label exchange"
sleep 60
python3 - > $S/count_fcc32_label_exchange_M2.txt 2>&1 <<'PY'
import platform, time, numba, necklace_enum as te
print(f"# {platform.platform()}, Python {platform.python_version()}, necklace-enum {te.__version__}, Numba {numba.__version__}")
P = te.fcc(2)
te.count_structures(P, 8, label_exchange=True)          # compile first
t = time.perf_counter()
c = te.count_structures(P, 32, label_exchange=True)[32]
print(f"binary fcc, n=32, up to label exchange: {c:,d} structures in {time.perf_counter() - t:.1f} s (one run, after compilation)")
PY

echo "== 5. enumlib, fcc 42:6 at n = 48"
sleep 60
W=$(mktemp -d)
cat > "$W/struct_enum.in" <<'IN'
fcc
bulk
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
 2 -nary case
 1
 0.0 0.0 0.0 0/1
 48 48
0.10000000E-06
full list of labelings
# Concentration restrictions
42 42 48
6 6 48
IN
ln -s /dev/null "$W/struct_enum.out"
{
  echo "# enumlib $(git -C "$ENUMLIB" log -1 --format='%h' ), fcc 42:6 at n = 48, output to /dev/null, stopped after 3600 s"
  cat "$W/struct_enum.in"
  echo "# ---- enum.x stdout and stderr, then /usr/bin/time -l"
  (cd "$W" && perl -e 'alarm shift; exec @ARGV' 3600 /usr/bin/time -l "$ENUMX" 2>&1)
  echo "# exit status: $?"
} > $S/enumlib_fcc42-6_n48_M2.txt 2>&1

echo "finished $(date)"
