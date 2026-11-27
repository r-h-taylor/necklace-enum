# Superseded benchmark logs

These timing runs were affected by an orphaned `enum.x` process (left over from an earlier memory-measurement
timeout) that ran in the background for many hours and slowed every code by 15-25%. They are kept for the record
and are not used in the paper. The process was found and killed, the timeout handling was fixed (the benchmark scripts, not included in this repository, were changed to kill the whole process group), and all timings were rerun on an
idle, externally air-cooled machine: see `../*_timings_M2.txt`.

Also kept here for the record, not used in the paper:

- `memory_M2.txt`: an early memory run with an earlier necklace-enum version (higher interpreter baseline and
  listing memory). Superseded by `../memory_big_M2.txt`, `../memory_dsenum_enumlibjl_M2.txt` and
  `../memory_enumlibjl050_M2.txt`.
- `enumlibjl_timings_M2.txt`: Enumlib.jl 0.4.2. Superseded by `../enumlibjl_v050_timings_M2.txt` (0.5.0).
