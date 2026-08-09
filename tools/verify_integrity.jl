# Regression checks for the persistent-cache and production-resume contracts.
# These tests only write below temporary directories and never touch real caches.

include(joinpath(@__DIR__, "..", "src", "gen_production.jl"))

function expect_rejected(f)
    rejected = false
    try
        f()
    catch
        rejected = true
    end
    @assert rejected
end

function verify_cache_integrity()
    mktempdir() do dir
        cd(dir) do
            empty!(_cache)
            key = ("integrity", 7)
            obj = Dict("values" => [1.0, 2.0], "ok" => true)
            cache_put(key, obj)
            fname = cache_file(key)
            @assert isfile(fname)
            @assert occursin(CACHE_SOURCE_FINGERPRINT, fname)

            empty!(_cache)
            @assert disk_cached(() -> error("unexpected rebuild"), key) == obj

            envelope = deserialize(fname)
            payload = copy(envelope.payload)
            payload[1] = xor(payload[1], 0x01)
            damaged = merge(envelope, (payload=payload,))
            expect_rejected(() -> cache_unwrap(damaged, key))

            # The public read path must rebuild a corrupt file, not merely throw.
            serialize(fname, damaged)
            empty!(_cache)
            rebuilt = Ref(false)
            got = disk_cached(() -> (rebuilt[] = true; Dict("rebuilt" => true)), key)
            @assert rebuilt[] && got["rebuilt"]
        end
    end
end

function sample_partial_row()
    F = zeros(length(S_GRID))
    F[1] = 1.0
    return Dict{String,Any}(
        "e0_keV" => 200.0, "u" => 2.0, "F" => F, "N0" => 1.0,
        "sigma_own_nm2" => 1.0, "sigma_bote_nm2" => 1.0,
        "s_cert_A_inv" => S_GRID[1],
        "tail" => Dict{String,Any}(
            "kind" => TAIL_KIND_BOUND, "source" => "integrity-test",
            "eps" => 1e-6, "valid_to" => S_GRID[1]),
        "diag" => Dict{String,Any}(
            "mres" => 1e-8, "badL" => 0, "rtail" => 1e-8,
            "ortho_c" => 1e-8, "retried" => 0))
end

function verify_checkpoint_integrity()
    row = sample_partial_row()
    context = production_context_sha256(QUICK_SETTINGS, PRESC_V4)
    record = checkpoint_record(row, context)
    @assert checkpoint_row(record, context) === row

    bad_value = deepcopy(record)
    bad_value["row"]["N0"] = 2.0
    expect_rejected(() -> checkpoint_row(bad_value, context))

    bad_context = deepcopy(record)
    bad_context["context_sha256"] = "wrong"
    expect_rejected(() -> checkpoint_row(bad_context, context))

    bad_diag = deepcopy(record)
    bad_diag["row"]["diag"]["badL"] = 0.5
    bad_diag["row_sha256"] = checkpoint_sha256(bad_diag["row"])
    expect_rejected(() -> checkpoint_row(bad_diag, context))

    expect_rejected(() -> checkpoint_row(row, context)) # old direct-row format

    mktempdir() do dir
        append_partial(dir, "K", 6, row, context)
        loaded = load_partial(dir, "K", 6, context)
        @assert length(loaded) == 1
        @assert loaded[200.0]["F"] == row["F"]
        @assert loaded[200.0]["diag"]["badL"] isa Int
        @assert isempty(load_partial(dir, "K", 6, "different-context"))
    end
end

verify_cache_integrity()
verify_checkpoint_integrity()
println("cache/checkpoint integrity: PASS")
