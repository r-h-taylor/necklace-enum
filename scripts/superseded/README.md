# Superseded benchmark logs

These timing runs were affected by an orphaned `enum.x` process (left over from an earlier memory-measurement
timeout) that ran in the background for many hours and slowed every code by 15-25%. They are kept for the record
and are not used in the paper. The process was found and killed, the timeout handling was fixed (the benchmark
scripts were changed to kill the whole process group), and all timings were rerun on an
idle, externally air-cooled machine: see `../*_timings_M2.txt`.
