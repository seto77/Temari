# cache_race_test.jl — persistent-cache publication/repair race regression
#
# This is intentionally a multi-process test.  A thread-only test cannot prove
# the filesystem create-if-absent property that protects a shared cwd.
# Everything is created below mktempdir(); the repository atom_cache is never
# read, removed, or written.
#
#   julia +1.11 -t 4 --gcthreads=1 --startup-file=no tools/cache_race_test.jl

using Serialization

include(joinpath(@__DIR__, "..", "src", "ionization.jl"))

const RACE_KEY = ("d", 999, 1, 0, -1, "race-test")
const LEGACY_KEY = ("cache-race-legacy", 1)

marker(path::String, text::String="ok") = open(path, "w") do io
    write(io, text)
end

function wait_until(pred, what::String; timeout_s::Float64=60.0)
    deadline = time() + timeout_s
    while !pred()
        time() < deadline || error("timeout waiting for $what")
        sleep(0.02)
    end
end

function child_cmd(args::Vector{String})
    Cmd(vcat(Base.julia_cmd().exec,
             ["-t", "1", "--gcthreads=1", @__FILE__], args))
end

function wait_success(processes)
    for p in processes
        wait(p)
        success(p) || error("child process failed: $p")
    end
end

function race_candidate(id::Int)
    # A recognized bound-state key makes every losing child run the semantic
    # validator, not only the envelope/checksum checks.  The payload is large
    # enough for independent serialization while still cheap next to real SCF.
    r = collect(range(1.0e-7, 1.0; length=8192))
    return (-Float64(id), r, fill(Float64(id), length(r)), 0.01)
end

function run_first_wins_child(work::String, id::Int)
    cd(work)
    marker(joinpath(work, "ready_$id"))
    wait_until(() -> isfile(joinpath(work, "start")), "race start")
    got = cache_put(RACE_KEY, race_candidate(id))
    serialize(joinpath(work, "returned_$id.jls"), got)
end

function run_sequential_child(work::String, id::Int)
    cd(work)
    marker(joinpath(work, "seq_ready_$id"))
    if id == 1
        wait_until(() -> isfile(joinpath(work, "seq_start")), "sequential start")
    else
        wait_until(() -> isfile(joinpath(work, "seq_done_1")), "writer A completion")
    end
    got = cache_put(RACE_KEY, race_candidate(id))
    serialize(joinpath(work, "seq_returned_$id.jls"), got)
    marker(joinpath(work, "seq_done_$id"))
end

function run_stale_lock_child(work::String)
    cd(work)
    key = ("cache-race-stale", 1)
    fname = cache_file(key)
    mkpath(dirname(fname))
    lockdir = cache_repair_lock(fname)
    cache_try_repair_lock(lockdir) === :acquired || error("child did not acquire lock")
    marker(joinpath(work, "stale_ready"))
    sleep(3600) # parent terminates us without finally/release
end

function legacy_candidate(id::String)
    Dict{String,Any}("writer" => id, "payload" => fill(codeunit(id, 1), 1 << 18))
end

"The old force-overwrite primitive, retained only to prove the negative test."
function legacy_publish(key::Tuple, obj)
    fname = cache_file(key)
    mkpath(dirname(fname))
    tmp = fname * ".legacy.$(getpid()).$(time_ns())"
    serialize(tmp, cache_envelope(key, obj))
    mv(tmp, fname; force=true)
    return cache_read_valid(fname, key)
end

function run_legacy_child(work::String, id::String)
    cd(work)
    marker(joinpath(work, "ready_$id"))
    if id == "A"
        wait_until(() -> isfile(joinpath(work, "start")), "legacy start")
        got = legacy_publish(LEGACY_KEY, legacy_candidate("A"))
        got["writer"] == "A" || error("first legacy writer did not publish A")
        marker(joinpath(work, "first_done"))
    else
        wait_until(() -> isfile(joinpath(work, "first_done")), "first legacy writer")
        first = cache_read_valid(cache_file(LEGACY_KEY), LEGACY_KEY)
        first["writer"] == "A" || error("second writer did not observe A")
        got = legacy_publish(LEGACY_KEY, legacy_candidate("B"))
        got["writer"] == "B" || error("second legacy writer did not overwrite A")
        marker(joinpath(work, "second_done"))
    end
end

function test_legacy_negative(root::String)
    work = joinpath(root, "legacy-negative")
    mkpath(work)
    processes = [run(child_cmd(["--legacy-child", work, id]); wait=false)
                 for id in ("A", "B")]
    wait_until(() -> all(isfile(joinpath(work, "ready_$id")) for id in ("A", "B")),
               "legacy children")
    marker(joinpath(work, "start"))
    wait_success(processes)
    cd(work) do
        final = cache_read_valid(cache_file(LEGACY_KEY), LEGACY_KEY)
        @assert final["writer"] == "B"
    end
    println("  PASS negative: old mv(force=true) let process B overwrite process A")
end

function test_atomic_first_wins(root::String)
    work = joinpath(root, "atomic-first-wins")
    mkpath(work)
    nchild = 8
    processes = [run(child_cmd(["--first-wins-child", work, string(id)]); wait=false)
                 for id in 1:nchild]
    wait_until(() -> all(isfile(joinpath(work, "ready_$id")) for id in 1:nchild),
               "first-wins children")
    marker(joinpath(work, "start"))
    wait_success(processes)
    cd(work) do
        empty!(_cache)
        winner = cache_read_valid(cache_file(RACE_KEY), RACE_KEY)
        returned = [deserialize(joinpath(work, "returned_$id.jls")) for id in 1:nchild]
        @assert all(x -> x == winner, returned)
        winner_id = round(Int, -winner[1])
        @assert count(id -> id != winner_id, 1:nchild) == nchild - 1
        leftovers = filter(name -> occursin(".tmp.", name) || endswith(name, ".repair.lock"),
                           readdir(CACHE_DIR))
        @assert isempty(leftovers)
        println("  PASS regression: writer $winner_id won; " *
                "all $(nchild - 1) losers returned the deserialized/validated winner")
    end
end

function test_sequential_first_wins(root::String)
    work = joinpath(root, "sequential-first-wins")
    mkpath(work)
    processes = [run(child_cmd(["--sequential-child", work, string(id)]); wait=false)
                 for id in 1:2]
    wait_until(() -> all(isfile(joinpath(work, "seq_ready_$id")) for id in 1:2),
               "sequential children")
    marker(joinpath(work, "seq_start"))
    wait_success(processes)
    cd(work) do
        empty!(_cache)
        disk = cache_read_valid(cache_file(RACE_KEY), RACE_KEY)
        a_return = deserialize(joinpath(work, "seq_returned_1.jls"))
        b_return = deserialize(joinpath(work, "seq_returned_2.jls"))
        @assert disk == race_candidate(1)
        @assert a_return == disk
        @assert b_return == disk
    end
    println("  PASS sequential: A completed first; later B returned A and disk remained A")
end

function test_stale_lock_recovery(root::String)
    work = joinpath(root, "stale-lock")
    mkpath(work)
    p = run(child_cmd(["--stale-lock-child", work]); wait=false)
    wait_until(() -> isfile(joinpath(work, "stale_ready")), "stale-lock child")
    kill(p)
    wait(p)
    cd(work) do
        key = ("cache-race-stale", 1)
        lockdir = cache_repair_lock(cache_file(key))
        @assert isdir(lockdir)
        owner = cache_parse_lock_record(cache_repair_owner(lockdir))
        @assert owner !== nothing
        @assert cache_process_state(owner.pid) === :dead
        # Production requires five minutes; zero is an explicit test hook that
        # still requires the recorded owner PID to be dead.
        @assert cache_try_reap_stale_lock(lockdir; stale_s=0.0)
        @assert !ispath(lockdir)
        @assert !ispath(cache_reap_claim(lockdir))
        got = cache_put(key, Dict{String,Any}("after_crash" => true))
        @assert got["after_crash"] === true
        @assert cache_read_valid(cache_file(key), key) == got
    end
    println("  PASS stale-lock: killed owner was proven dead and its old lock reclaimed")
end

function quarantine_files()
    isdir(CACHE_DIR) || return String[]
    filter(path -> occursin(".corrupt.", basename(path)), readdir(CACHE_DIR; join=true))
end

function test_corrupt_quarantine(root::String)
    work = joinpath(root, "corrupt-repair")
    mkpath(work)
    cd(work) do
        empty!(_cache)
        key = ("cache-race-corrupt", 1)
        fname = cache_file(key)
        mkpath(dirname(fname))
        broken = UInt8[0x54, 0x65, 0x6d, 0x61, 0x72, 0x69, 0x00, 0xff]
        write(fname, broken)
        candidate = Dict{String,Any}("writer" => "repair", "ok" => true)
        got = cache_put(key, candidate)
        @assert got == candidate
        @assert cache_read_valid(fname, key) == candidate
        qs = quarantine_files()
        @assert length(qs) == 1
        @assert read(only(qs)) == broken
        @assert !isdir(cache_repair_lock(fname))
    end
    println("  PASS repair: corrupt bytes were preserved exactly and replaced under lock")
end

function test_semantic_quarantine(root::String)
    work = joinpath(root, "semantic-repair")
    mkpath(work)
    cd(work) do
        empty!(_cache)
        cfg = NumericsConfig()
        key = neutral_cache_key(2, false, X_ALPHA, :xalpha, cfg)
        fname = cache_file(key)
        mkpath(dirname(fname))
        # Envelope/key/checksum are internally consistent, but the payload is
        # semantically impossible for a neutral SCF key.
        serialize(fname, cache_envelope(key, Dict("not" => "an SCFAtom")))
        bad_bytes = read(fname)
        try
            cache_read_valid(fname, key)
            error("semantic corruption was accepted")
        catch err
            occursin("SCF cache object type mismatch", sprint(showerror, err)) || rethrow()
        end
        good = build_neutral(2)
        got = cache_put(key, good)
        @assert got isa SCFAtom && got.z == 2
        @assert cache_read_valid(fname, key).z == 2
        qs = quarantine_files()
        @assert length(qs) == 1
        @assert read(only(qs)) == bad_bytes
        @assert !isdir(cache_repair_lock(fname))
    end
    println("  PASS semantics: valid-checksum/wrong-object payload was rejected and preserved")
end

function test_uncertain_path_fails_closed(root::String)
    work = joinpath(root, "uncertain-path")
    mkpath(work)
    cd(work) do
        empty!(_cache)
        key = ("cache-race-uncertain", 1)
        fname = cache_file(key)
        mkpath(fname) # unexpected object at the exact public cache name
        sentinel = joinpath(fname, "keep.txt")
        marker(sentinel, "do-not-change")
        candidate = Dict{String,Any}("memory" => true)
        got = cache_put(key, candidate)
        @assert got == candidate
        @assert isdir(fname)
        @assert read(sentinel, String) == "do-not-change"
        @assert !isdir(cache_repair_lock(fname))
    end
    println("  PASS fail-closed: uncertain destination was untouched; memory candidate returned")
end

function main()
    mktempdir() do root
        println("isolated root: ", root)
        test_legacy_negative(root)
        test_sequential_first_wins(root)
        test_atomic_first_wins(root)
        test_stale_lock_recovery(root)
        test_corrupt_quarantine(root)
        test_semantic_quarantine(root)
        test_uncertain_path_fails_closed(root)
    end
    println("cache race/quarantine: ALL PASS")
end

if !isempty(ARGS) && ARGS[1] == "--first-wins-child"
    run_first_wins_child(ARGS[2], parse(Int, ARGS[3]))
elseif !isempty(ARGS) && ARGS[1] == "--sequential-child"
    run_sequential_child(ARGS[2], parse(Int, ARGS[3]))
elseif !isempty(ARGS) && ARGS[1] == "--stale-lock-child"
    run_stale_lock_child(ARGS[2])
elseif !isempty(ARGS) && ARGS[1] == "--legacy-child"
    run_legacy_child(ARGS[2], ARGS[3])
else
    main()
end
