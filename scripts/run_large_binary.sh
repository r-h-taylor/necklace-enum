#!/bin/sh
# Binary fcc over all compositions at the largest sizes, with a one-hour limit per run:
#   1. the present method listing as integer arrays at n = 30 and 32 (time and peak memory);
#   2. dsenum counting on the larger problems of Table 2.
# Same protocol as the paper: one core, otherwise idle, externally cooled, 60 s pause before
# each timed run.
#
#   cd ~/Projects/necklace-enum
#   nohup caffeinate -i sh scripts/run_large_binary.sh < /dev/null > /dev/null 2>&1 &
#
# Progress: cat scripts/sweep_binary_fcc_3032_M2.txt scripts/memory_big_dsenum_1h_M2.txt
cd "$(dirname "$0")/.." || exit 1

python3 scripts/sweep_binary_fcc.py --codes necklace --sizes 30,32 --extra "" \
    --cooldown 60 --timeout 3600 > scripts/sweep_binary_fcc_3032_M2.txt 2>&1

sleep 60
python3 scripts/memory_compare.py --codes dsenum --only-big --timeout 3600 \
    > scripts/memory_big_dsenum_1h_M2.txt 2>&1
