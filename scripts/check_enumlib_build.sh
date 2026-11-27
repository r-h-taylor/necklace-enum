#!/bin/bash
# Check how enum.x was built (compiler, architecture, optimisation, debug flags) and,
# with --memory, re-measure enumlib's peak memory on the fcc 24:8, n=32 problem.
#
#   bash scripts/check_enumlib_build.sh ~/Projects/enumlib            # build checks only
#   bash scripts/check_enumlib_build.sh ~/Projects/enumlib --memory   # also the memory run (~5 min)

EL="${1:-$HOME/Projects/enumlib}"
EL="${EL/#\~/$HOME}"
X="$EL/src/enum.x"
[ -x "$X" ] || { echo "no executable at $X"; exit 1; }

echo "== enum.x"
ls -l "$X"
file "$X"                                   # should say arm64

echo; echo "== Fortran runtime it links to (tells which gfortran built it)"
otool -L "$X" | grep -i -E "gfortran|quadmath|gcc" || echo "(libgfortran linked statically or not found)"
LIB=$(otool -L "$X" | awk '/libgfortran/{print $1; exit}')
if [ -n "$LIB" ]; then
  # Homebrew paths look like /opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (gcc formula)
  for G in /opt/homebrew/bin/gfortran /usr/local/bin/gfortran; do
    [ -x "$G" ] && echo "$G: $($G --version | head -1)  [lib dir: $($G -print-file-name=libgfortran.dylib)]"
  done
fi
echo "gfortran first on PATH now: $(command -v gfortran) -> $(gfortran --version 2>/dev/null | head -1)"

echo; echo "== Compiler recorded in the object files (if present)"
for o in "$EL"/src/derivative_structure_generator.o "$EL"/symlib/src/symmetry.o; do
  [ -f "$o" ] || continue
  v=$(strings -a "$o" | grep -m1 -E "GNU Fortran|GCC: \(")
  echo "$(basename "$o"): ${v:-no version string (normal on macOS)}"
done

echo; echo "== Optimisation and debug flags"
echo "src/Makefile flags for gfortran:"; grep -n -E "^\s*FFLAGS|O3|DEBUG" "$EL/src/Makefile" | head -8
for a in "$EL"/symlib/src/libutils.a "$EL"/symlib/src/libsym.a; do
  [ -f "$a" ] || continue
  dbg=$(otool -l "$a" 2>/dev/null | grep -c __debug_info)
  bnd=$(nm "$a" 2>/dev/null | grep -c runtime_error_at)
  echo "$(basename "$a"): debug-info sections: $dbg  bounds-check calls: $bnd   (0 and 0 = optimised build; >0 = -g / -fbounds-check)"
done
dbg=$(otool -l "$X" | grep -c __DWARF); bnd=$(nm "$X" 2>/dev/null | grep -c runtime_error_at)
echo "enum.x: DWARF segments: $dbg  bounds-check calls: $bnd"

if [ "$2" = "--memory" ]; then
  echo; echo "== Peak memory, fcc 24:8 at n=32 (output to /dev/null)"
  W=$(mktemp -d)
  cat > "$W/struct_enum.in" <<'EOF'
fcc
bulk
0.5 0.5 0.0
0.5 0.0 0.5
0.0 0.5 0.5
 2 -nary case
 1
 0.0 0.0 0.0 0/1
 32 32
0.10000000E-06
full list of labelings
# Concentration restrictions
24 24 32
8 8 32
EOF
  ln -s /dev/null "$W/struct_enum.out"
  (cd "$W" && /usr/bin/time -l "$X" > /dev/null 2> time.log)
  grep -E "real|maximum resident set size|peak memory footprint" "$W/time.log"
  awk '/maximum resident set size/{printf "peak RSS = %.2f GB\n", $1/1e9}' "$W/time.log"
fi
