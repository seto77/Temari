# jobq_rows_sigma.jl の LPT / sidecar / tag 較正 gate の負・回帰テスト (260823Cl)

include(joinpath(@__DIR__, "..", "..", "jobq_rows_sigma.jl"))

const TEST_CODE_SHA = "a"^64
const TEST_CERT_FP = "b"^16

function run_tests()
    passed = Ref(0); failed = Ref(0)
    function ok(cond, what)
        if cond
            passed[] += 1; println("  ok   ", what)
        else
            failed[] += 1; println("  FAIL ", what)
        end
    end
    function fails(f, needle, what)
        hit = false; msg = ""
        try
            f()
        catch e
            msg = sprint(showerror, e); hit = occursin(needle, msg)
        end
        ok(hit, what * (hit ? "" : " (message=$msg)"))
    end
    row_of_arg(a) = (Int(a["rows"][1][1]), String(a["rows"][1][2]), Float64(a["rows"][1][3]))

    mktempdir() do tmp
        gate_args = joinpath(tmp, "gate.args.json"); gate_est = joinpath(tmp, "gate.est.json")
        common = ["--profile", "deep", "--rule", "v4", "--code-sha256", TEST_CODE_SHA, "--cert-fp", TEST_CERT_FP]

        println("[1] sentinel gate: 11 行・生代理値 LPT・裸の args")
        ok(main_rows(vcat(common, ["--wave", "sentinel", "--out", gate_args, "--est-out", gate_est])) == 0,
           "sentinel generator exit 0")
        ga = parse_json_file(gate_args); ge = parse_json_file(gate_est)
        ok(length(ga) == 11 && Int(ge["n_jobs"]) == 11 && Int(ge["n_sentinel"]) == 11,
           "canonical sentinel 11 行 = 11 票")
        ok(all(sort(collect(keys(a))) == ["rows", "rule"] && length(a["rows"]) == 1 for a in ga),
           "args は rule+rows だけ・1 行 1 票 (est_min を混ぜない)")
        ok(!occursin("est_min", read(gate_args, String)) && haskey(ge["jobs"][1], "est_min"),
           "est_min は別 sidecar だけ")
        ok(row_of_arg(ga[1]) == (20, "M1", 400.0), "生代理 LPT の先頭は Ca M1@400")
        raw_gate = Float64[j["raw_cost"] for j in ge["jobs"]]
        ok(issorted(raw_gate; rev=true), "sentinel は raw_cost の降順")
        ok(ge["code_sha256"] == TEST_CODE_SHA && ge["cert_fp"] == TEST_CERT_FP,
           "gate sidecar に code_sha256 / cert_fp")

        println("[2] fail-closed: wave / 粒度 / calibration")
        fails(() -> main_rows(vcat(common, ["--out", joinpath(tmp, "x"), "--est-out", joinpath(tmp, "xe")])),
              "--wave", "deep で wave 省略を拒否")
        fails(() -> main_rows(vcat(common, ["--wave", "sentinel", "--group", "channel",
                                            "--out", joinpath(tmp, "x2"), "--est-out", joinpath(tmp, "xe2")])),
              "1 行 1 票", "deep の channel 束ねを拒否")
        fails(() -> main_rows(["--profile", "deep", "--wave", "sentinel", "--rule", "v4",
                               "--tags", "K", "--code-sha256", TEST_CODE_SHA, "--cert-fp", TEST_CERT_FP,
                               "--out", joinpath(tmp, "subset.args"), "--est-out", joinpath(tmp, "subset.est")]),
              "canonical 9 tags", "deep を tag 部分集合へ分割せず最終 1,583 行を守る")
        full_args = joinpath(tmp, "full.args.json"); full_est = joinpath(tmp, "full.est.json")
        fails(() -> main_rows(vcat(common, ["--wave", "full", "--out", full_args, "--est-out", full_est])),
              "--calibration", "full は calibration 無しを拒否")
        ok(!isfile(full_args) && !isfile(full_est), "拒否時に args / sidecar を残さない")

        pilot_args = joinpath(tmp, "pilot.args.json")
        ok(main_rows(["--profile", "pilot", "--rule", "v4", "--out", pilot_args]) == 0 &&
           !isempty(parse_json_file(pilot_args)) && !occursin("est_min", read(pilot_args, String)),
           "既存 pilot CLI は sidecar 無しの軽量生成として後方互換")

        println("[3] gate result + manifest から tag 別較正")
        result_dir = joinpath(tmp, "results"); mkpath(result_dir)
        # M5 を重めにするなど tag 係数を変え、full の cross-tag 順序が本当に較正されることも試す。
        factors = Dict("K" => 1.00, "L1" => 1.10, "M1" => 1.20, "M2" => 1.30,
                       "M3" => 1.40, "M4" => 1.50, "M5" => 1.60)
        for j in ge["jobs"]
            jobseq = Int(j["jobseq"]); z = Int(j["z"]); tag = String(j["tag"]); e0 = Float64(j["e0_keV"])
            raw = Float64(j["raw_cost"]); total = raw * factors[tag] * 60
            outname = @sprintf("temari_sigma_gate_test_lane%06d001.jsonl", jobseq)
            path = joinpath(result_dir, outname)
            ws = window_list(Float64(j["eps_max_eV"])); nw = length(ws)
            open(path, "w") do io
                jobseq == 1 && println(io, "{\"z\":$z,\"tag\":", _jq_str(tag),
                                       ",\"e0_keV\":", repr(e0), ",\"error\":\"transient fixture\",\"cert_fp\":",
                                       _jq_str(TEST_CERT_FP), "}")
                for (iw, w) in enumerate(ws)
                    # 実ファイルと同じく 1 JSON object / 1 line。較正に必要な欄だけの fixture。
                    rec = Dict{String,Any}(
                        "z" => z, "tag" => tag, "e0_keV" => e0, "window_id" => w[1],
                        "n_windows_in_row" => nw, "cert_fp" => TEST_CERT_FP,
                        "row_elapsed_s" => total * iw / nw, "out_of_domain" => !w[4])
                    if w[4]
                        merge!(rec, Dict{String,Any}(
                            "rule" => _expected_rule_name("v4"), "rule_version" => "v4",
                            "rule_config" => string(rule_config(_rule_object("v4"))),
                            "oracle" => _expected_oracle("v4"), "pass" => true))
                    end
                    println(io, jsonl_line_v2(rec))
                end
            end
            manifest = Dict{String,Any}(
                "schema" => 1, "campaign" => "temari_sigma_gate_test", "jobseq" => jobseq,
                "claim_epoch" => 1, "task" => "temari.certify_sigma_v2", "code_sha256" => TEST_CODE_SHA,
                "hostname" => "host-a",
                "outname" => outname, "result_sha256" => bytes2hex(sha256(read(path))),
                "task_info" => Dict("cert_fp" => [TEST_CERT_FP], "rows_done" => 1,
                                    "error_lines" => (jobseq == 1 ? 1 : 0)))
            open(path * ".manifest.json", "w") do io; write_json(io, manifest); write(io, '\n'); end
        end

        println("[3a] gate result の全仕様内窓に規則・oracle・domain が必要")
        sample = joinpath(result_dir, sort([f for f in readdir(result_dir) if endswith(f, ".jsonl")])[1])
        sm = parse_json_file(sample * ".manifest.json"); sj = ge["jobs"][Int(sm["jobseq"])]
        sample_row = (Int(sj["z"]), String(sj["tag"]), Float64(sj["e0_keV"]))
        function mutate_sample(path, mutate)
            rows = [_parse_json_line(line, sample, i) for (i, line) in enumerate(eachline(sample))]
            idx = findfirst(d -> !haskey(d, "error") && get(d, "out_of_domain", false) === false, rows)
            idx === nothing && error("fixture に仕様内窓が無い")
            mutate(rows[idx])
            open(path, "w") do io
                for d in rows; println(io, jsonl_line_v2(d)); end
            end
        end
        missing_oracle = joinpath(tmp, "missing-oracle.jsonl")
        mutate_sample(missing_oracle, d -> delete!(d, "oracle"))
        fails(() -> _complete_attempt(missing_oracle, sample_row, TEST_CERT_FP, "v4"),
              "揃っていない", "1 窓だけ oracle が無い gate 結果を拒否")
        false_ood = joinpath(tmp, "false-out-of-domain.jsonl")
        mutate_sample(false_ood, d -> (d["out_of_domain"] = true))
        fails(() -> _complete_attempt(false_ood, sample_row, TEST_CERT_FP, "v4"),
              "canonical domain", "仕様内窓を out_of_domain と偽る結果を拒否")

        cal = joinpath(tmp, "calibration.json")
        fails(() -> main_calibration(["--sentinel-est", gate_est, "--results-dir", result_dir, "--out", cal]),
              "tag L2", "sentinel に無い L2/L3 は明示 fallback 無しでは拒否")
        ok(!isfile(cal), "較正失敗時に sidecar を残さない")
        fails(() -> main_calibration(["--sentinel-est", gate_est, "--results-dir", result_dir,
                                      "--fallback", "K=L1", "--fallback", "L2=L1", "--fallback", "L3=L1",
                                      "--out", joinpath(tmp, "measured-fallback.json")]),
              "実測済み tag K", "実測 tag を fallback で黙って上書きしない")
        ok(main_calibration(["--sentinel-est", gate_est, "--results-dir", result_dir,
                             "--fallback", "L2=L1", "--fallback", "L3=L1", "--out", cal]) == 0,
           "明示 fallback つき calibration exit 0")
        c = parse_json_file(cal); tc = c["tag_calibration"]
        ok(all(isapprox(Float64(tc[t]["minutes_per_proxy"]), factors[t]; rtol=2e-15) for t in keys(factors)),
           "実測 elapsed_min/raw_cost の tag 中央値")
        ok(tc["L2"]["source"] == "fallback:L1" && tc["L3"]["source"] == "fallback:L1",
           "未観測 tag の由来を sidecar に明記")
        ok(c["source_campaign"] == "temari_sigma_gate_test" && c["code_sha256"] == TEST_CODE_SHA &&
           c["cert_fp"] == TEST_CERT_FP && length(c["observations"]) == 11,
           "較正 sidecar に campaign / code / cert / 11 観測の来歴")
        ok(sum(Int(o["recovered_error_lines"]) for o in c["observations"]) == 1,
           "完結 attempt 前の一時 error は queue verifier と同じ規則で許し来歴へ残す")
        ok(c["hosts"] == ["host-a"] && Float64(c["host_slowdown"]["host-a"]) == 1.0 &&
           all(Float64(o["elapsed_s"]) == Float64(o["normalized_elapsed_s"]) for o in c["observations"]),
           "単一 hostname は slowdown=1 として明記")

        println("[3b] full は calibration の埋込み来歴を再検算")
        function write_obj(path, obj)
            open(path, "w") do io; write_json(io, obj); write(io, '\n'); end
        end
        function reject_cal(obj, stem, needle, what)
            path = joinpath(tmp, stem * ".json"); a = joinpath(tmp, stem * ".args.json"); e = joinpath(tmp, stem * ".est.json")
            write_obj(path, obj)
            fails(() -> main_rows(vcat(common, ["--wave", "full", "--calibration", path,
                                                "--out", a, "--est-out", e])), needle, what)
            ok(!isfile(a) && !isfile(e), what * " (成果物を残さない)")
        end
        no_observations = parse_json_file(cal); delete!(no_observations, "observations")
        reject_cal(no_observations, "cal-no-observations", "observations", "observation 無しの捏造 calibration を拒否")
        bad_measured = parse_json_file(cal)
        bad_measured["tag_calibration"]["K"]["minutes_per_proxy"] =
            Float64(bad_measured["tag_calibration"]["K"]["minutes_per_proxy"]) * 1.01
        reject_cal(bad_measured, "cal-bad-measured", "measured coefficient", "実測 tag 係数の改変を拒否")
        bad_fallback = parse_json_file(cal)
        bad_fallback["tag_calibration"]["L2"]["minutes_per_proxy"] =
            Float64(bad_fallback["tag_calibration"]["L2"]["minutes_per_proxy"]) * 1.01
        reject_cal(bad_fallback, "cal-bad-fallback", "fallback coefficient", "fallback 係数の改変を拒否")
        missing_result = parse_json_file(cal); pop!(missing_result["result_files"])
        reject_cal(missing_result, "cal-missing-result", "件数", "result/manifest 来歴の欠落を拒否")

        println("[3c] 複数 hostname の機速正規化は fail-closed")
        first_result = joinpath(result_dir, sort([f for f in readdir(result_dir) if endswith(f, ".jsonl")])[1])
        first_manifest = first_result * ".manifest.json"
        mm = parse_json_file(first_manifest); mm["hostname"] = "host-b"
        open(first_manifest, "w") do io; write_json(io, mm); write(io, '\n'); end
        multi_base = ["--sentinel-est", gate_est, "--results-dir", result_dir,
                      "--fallback", "L2=L1", "--fallback", "L3=L1"]
        fails(() -> main_calibration(vcat(multi_base, ["--out", joinpath(tmp, "multi-missing.json")])),
              "全 host", "複数 hostname を slowdown 無しで平均しない")
        fails(() -> main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=1", "--out", joinpath(tmp, "multi-one.json")])),
              "host-b", "複数 hostname の slowdown 欠落を列挙")
        fails(() -> main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=0", "--out", joinpath(tmp, "multi-zero.json")])),
              "正の有限値", "host slowdown 0 を拒否")
        fails(() -> main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=NaN", "--out", joinpath(tmp, "multi-nan.json")])),
              "正の有限値", "host slowdown NaN を拒否")
        fails(() -> main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=1", "--host-slowdown", "HOST-A=2",
                                                    "--out", joinpath(tmp, "multi-dup.json")])),
              "重複", "host slowdown の大文字小文字を跨ぐ重複を拒否")
        fails(() -> main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=1", "--host-slowdown", "host-b=2",
                                                    "--host-slowdown", "ghost=1", "--out", joinpath(tmp, "multi-unknown.json")])),
              "結果に無い", "結果に無い host slowdown を拒否")
        multi_cal = joinpath(tmp, "multi-valid.json")
        ok(main_calibration(vcat(multi_base, ["--host-slowdown", "host-a=1", "--host-slowdown", "host-b=2",
                                                   "--out", multi_cal])) == 0,
           "複数 hostname を全係数指定すれば較正")
        mc = parse_json_file(multi_cal); mbo = only(o for o in mc["observations"] if o["hostname"] == "host-b")
        ok(Float64(mbo["normalized_elapsed_s"]) == Float64(mbo["elapsed_s"]) / 2 &&
           occursin("elapsed_s/host_slowdown", mc["normalization_formula"]),
           "elapsed/slowdown で共通ホスト尺度へ正規化し式も記録")

        println("[4] full: 較正済 LPT・1,583 行・決定論・single final campaign 用")
        ok(main_rows(vcat(common, ["--wave", "full", "--calibration", cal,
                                   "--out", full_args, "--est-out", full_est])) == 0,
           "full generator exit 0")
        fa = parse_json_file(full_args); fe = parse_json_file(full_est)
        ok(length(fa) == 1583 && Int(fe["n_jobs"]) == 1583 && Int(fe["n_sentinel"]) == 11,
           "full = 1,583 行、sentinel 11 行も最終 campaign で再計算")
        ok(all(sort(collect(keys(a))) == ["rows", "rule"] && length(a["rows"]) == 1 for a in fa),
           "full args も rule+rows だけ・1 行 1 票")
        jobs = fe["jobs"]; estimates = Float64[j["est_min"] for j in jobs]
        ok(issorted(estimates; rev=true), "full は tag 較正後 est_min の降順")
        ok(all(isapprox(Float64(j["est_min"]), Float64(j["raw_cost"]) *
                        Float64(tc[String(j["tag"])]["minutes_per_proxy"]); rtol=2e-15) for j in jobs),
           "est_min = raw_cost × tag minutes/proxy")
        ok(all(row_of_arg(fa[i]) == (Int(jobs[i]["z"]), String(jobs[i]["tag"]), Float64(jobs[i]["e0_keV"]))
               for i in eachindex(fa)), "args jobseq と sidecar jobseq が一致")
        full_set = Set(_rowkey(row_of_arg(a)) for a in fa)
        canonical = JobqSigmaRow.(profile_rows("deep", copy(TAGS_V4), ""))
        ok(length(full_set) == 1583 && full_set == Set(_rowkey.(canonical)), "canonical deep 行集合を無欠落・無重複で保持")
        ca = findfirst(j -> Int(j["z"]) == 20 && j["tag"] == "M1" && Float64(j["e0_keV"]) == 400.0, jobs)
        ok(ca !== nothing && jobs[ca]["is_sentinel"] === true, "Ca M1@400 は full に残り sentinel と識別")
        raw_order, _ = lpt_rows(canonical, nothing)
        ok(_rowkey.(raw_order) != [_rowkey(row_of_arg(a)) for a in fa], "tag 係数が cross-tag 順序へ実際に効く")
        ok(row_of_arg(fa[1]) != canonical[1], "旧 tag 順 (先頭 K) の反 LPT へ退行していない")

        full_args2 = joinpath(tmp, "full2.args.json"); full_est2 = joinpath(tmp, "full2.est.json")
        ok(main_rows(vcat(common, ["--wave", "full", "--calibration", cal,
                                   "--out", full_args2, "--est-out", full_est2])) == 0 &&
           read(full_args) == read(full_args2) && read(full_est) == read(full_est2),
           "同じ入力から args / est sidecar が byte deterministic")

        println("[5] provenance の負テスト")
        tampered_gate = parse_json_file(gate_est)
        tampered_gate["jobs"][1]["raw_cost"] = Float64(tampered_gate["jobs"][1]["raw_cost"]) + 1.0
        tampered_gate_path = joinpath(tmp, "tampered-gate.est.json")
        open(tampered_gate_path, "w") do io; write_json(io, tampered_gate); write(io, '\n'); end
        fails(() -> main_calibration(["--sentinel-est", tampered_gate_path, "--results-dir", result_dir,
                                      "--fallback", "L2=L1", "--fallback", "L3=L1",
                                      "--out", joinpath(tmp, "tampered-gate-cal.json")]),
              "費用根拠", "gate sidecar の raw cost 改変を再計算で拒否")
        badcal = parse_json_file(cal); badcal["cert_fp"] = "c"^16
        badpath = joinpath(tmp, "bad-calibration.json")
        open(badpath, "w") do io; write_json(io, badcal); write(io, '\n'); end
        badargs = joinpath(tmp, "bad.args.json"); badest = joinpath(tmp, "bad.est.json")
        fails(() -> main_rows(vcat(common, ["--wave", "full", "--calibration", badpath,
                                            "--out", badargs, "--est-out", badest])),
              "cert_fp", "full は calibration の cert_fp 不一致を拒否")
        ok(!isfile(badargs) && !isfile(badest), "provenance 不一致で成果物を残さない")

        println("  order top3 = ", join(["$(jobs[i]["tag"]) Z$(Int(jobs[i]["z"])) @$(jobs[i]["e0_keV"])" for i in 1:3], ", "))
        println("  order last = ", "$(jobs[end]["tag"]) Z$(Int(jobs[end]["z"])) @$(jobs[end]["e0_keV"])")
        println("  Ca M1@400 final jobseq = ", ca)
    end
    println("jobq_rows_sigma_test: PASS $(passed[]) / FAIL $(failed[])")
    return failed[] == 0 ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(run_tests())
