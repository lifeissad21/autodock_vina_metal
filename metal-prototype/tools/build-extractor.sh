#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
VINA="$ROOT/../vendor/AutoDock-Vina/src/lib"

clang++ -std=c++17 -O2 \
  -I/opt/homebrew/include -I"$VINA" \
  "$ROOT/tools/extract_vina_model.cpp" \
  "$VINA/ad4cache.cpp" "$VINA/cache.cpp" "$VINA/conf_independent.cpp" \
  "$VINA/coords.cpp" "$VINA/grid.cpp" "$VINA/model.cpp" "$VINA/monte_carlo.cpp" \
  "$VINA/mutate.cpp" "$VINA/non_cache.cpp" "$VINA/parallel_mc.cpp" \
  "$VINA/parallel_progress.cpp" "$VINA/parse_pdbqt.cpp" "$VINA/quasi_newton.cpp" \
  "$VINA/quaternion.cpp" "$VINA/random.cpp" "$VINA/szv_grid.cpp" "$VINA/utils.cpp" \
  -L/opt/homebrew/lib -lboost_filesystem -lboost_thread -lboost_serialization \
  -lboost_program_options -o "$ROOT/tools/extract_vina_model"

