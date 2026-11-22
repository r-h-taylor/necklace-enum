# Time Enumlib.jl (Hart's Julia reimplementation of enumlib) on one problem.
#
#   julia scripts/time_enumlibjl.jl PARENT CELLS SPECIES COMPOSITION [--skip-resource-check]
#
#   PARENT       fcc | hcp
#   CELLS        number of primitive cells n (hcp: 2n sites)
#   SPECIES      number of atom types k
#   COMPOSITION  "all" for every composition, or atoms per cell, e.g. "24,8" or "10,5,5"
#
# Prints one line:  RESULT <structures> <seconds>   (or REFUSED / ERROR with a message).
# The time is for enumerate_structures alone, after a warm-up call on a small cell so that
# compilation is excluded. Called by scripts/compare_with_enumlibjl.py.

using Enumlib

parent_name = ARGS[1]
n = parse(Int, ARGS[2])
k = parse(Int, ARGS[3])
comp = ARGS[4] == "all" ? nothing : parse.(Int, split(ARGS[4], ","))
skip = "--skip-resource-check" in ARGS

if parent_name == "fcc"
    p = ParentLattice([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0])
    nsites = 1
elseif parent_name == "hcp"
    A = [1.0 -0.5 0.0; 0.0 sqrt(3) / 2 0.0; 0.0 0.0 sqrt(8 / 3)]
    p = ParentLattice(A, [[0.0, 0.0, 0.0], [1 / 3, 2 / 3, 1 / 2]])
    nsites = 2
else
    error("unknown parent $parent_name")
end
s = Sites(p, collect(0:k-1))

function run(nc, c)
    kw = c === nothing ? (;) : (; concentration = concentration_count(c; n_total = nc * nsites))
    enumerate_structures(p, s; supercells = VolumeRange(nc:nc), skip_resource_check = skip, kw...)
end

# warm-up on a small cell of the same kind (compilation excluded from the timing)
if comp === nothing
    run(2, nothing)
else
    run(k, nsites == 1 ? ones(Int, k) : 2 .* ones(Int, k))
end

try
    t = @elapsed e = run(n, comp)
    println("RESULT ", length(e), " ", t)
catch err
    msg = replace(sprint(showerror, err), '\n' => ' ')
    if occursin("memory", lowercase(msg)) || occursin("budget", lowercase(msg)) || occursin("cost", lowercase(msg))
        println("REFUSED ", first(msg, 200))
    else
        println("ERROR ", first(msg, 200))
    end
end
