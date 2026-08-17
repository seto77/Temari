# L0 JSON — 標準ライブラリに JSON が無いための最小パーサ / ライタ
#
# 用途は bote_salvat.json・reference_values.json の読み込みと --json 保存のみ。
# 数値ではないが、同じく「Python の標準品を自前で置き換えた小道具」なので L0 に置く。
#
# ⚠ 260818Cl 訂正: 上の「用途は … のみ」は成り立たない。出荷データセット
#   (prod_v5_jl / prod_factors_v1) の JSON 本体と本番チェックポイントもこの codec が
#   書き (gen_production.jl / gen_factors.jl の write_json)、公開 loader
#   (tools/factors_loader.jl) は parse_json_file に依存する。canonical は JSON なので
#   (docs/dataset_contract_2026-08-09.md)、ここは補助の小道具ではなく出荷物の
#   直列化器そのもの。下の 260809Cl エスケープ欠陥が重かったのもそのため。

# ---- 最小 JSON パーサ (自前の JSON 一式用。260818Cl: 「2 ファイル専用」ではない) ----
# 対応: object / array / 数値 / 文字列 / true / false / null。UTF-8 の
# マルチバイト文字 (日本語コメント等) を扱うためバイト列上で走査する。
# 標準ライブラリに JSON が無いための自前品。

const _WS = (0x20, 0x09, 0x0a, 0x0d)   # 空白, タブ, LF, CR

function _json_value(b::Vector{UInt8}, i::Int)
    while b[i] in _WS; i += 1; end
    c = b[i]
    if c == UInt8('{')
        obj = Dict{String,Any}()
        i += 1
        while true
            while b[i] in _WS; i += 1; end
            b[i] == UInt8('}') && return obj, i + 1
            key, i = _json_value(b, i)
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(':') || error("JSON: ':' expected at byte $i")
            val, i = _json_value(b, i + 1)
            obj[key] = val
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(',') && (i += 1)
        end
    elseif c == UInt8('[')
        arr = Any[]
        i += 1
        while true
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(']') && return arr, i + 1
            val, i = _json_value(b, i)
            push!(arr, val)
            while b[i] in _WS; i += 1; end
            b[i] == UInt8(',') && (i += 1)
        end
    elseif c == UInt8('"')
        buf = IOBuffer()
        i += 1
        while b[i] != UInt8('"')               # 0x22 は UTF-8 継続バイトと衝突しない
            if b[i] == UInt8('\\')
                # 260809Cl: \n と \t しか解いておらず、\r \b \f と \uXXXX が
                # **黙って文字そのもの**になっていた ("A" が "u0041" になる)。
                # writer 側 (`json_escape`) と往復が閉じるように揃えた
                i += 1
                ch = b[i]
                if ch == UInt8('n');     write(buf, 0x0a)
                elseif ch == UInt8('t'); write(buf, 0x09)
                elseif ch == UInt8('r'); write(buf, 0x0d)
                elseif ch == UInt8('b'); write(buf, 0x08)
                elseif ch == UInt8('f'); write(buf, 0x0c)
                elseif ch == UInt8('u')
                    cp = parse(UInt16, String(b[i+1:i+4]); base=16)
                    # ⚠ サロゲート対は扱わない (この codec が書く範囲では
                    #    制御文字しか \u で出さないため)。来たらエラーで止める
                    0xD800 <= cp <= 0xDFFF &&
                        error("JSON: サロゲート対は未対応 (byte $i)")
                    print(buf, Char(cp))
                    i += 4
                else                     # \" \\ \/ はその文字そのもの
                    write(buf, ch)
                end
            else
                write(buf, b[i])
            end
            i += 1
        end
        return String(take!(buf)), i + 1
    elseif c == UInt8('t')
        return true, i + 4
    elseif c == UInt8('f')
        return false, i + 5
    elseif c == UInt8('n')
        return nothing, i + 4
    else
        j = i
        while j <= length(b) && (UInt8('0') <= b[j] <= UInt8('9') ||
                                 b[j] in (UInt8('+'), UInt8('-'), UInt8('.'),
                                          UInt8('e'), UInt8('E')))
            j += 1
        end
        return parse(Float64, String(b[i:j-1])), j
    end
end

"JSON ファイルを読む (依存を持たないための最小実装。260818Cl: 出荷 JSON もこれで読む)"
parse_json_file(path::String) = _json_value(read(path), 1)[1]

# ---- 最小 JSON writer (--json 保存用) ----

"""JSON 文字列のエスケープ (260809Cl 追加)。

⚠ **それまで writer は一切エスケープしていなかった** — `"` を含む文字列を書くと
**自分の `parse_json_file` すら読み戻せない JSON** が出る (`tools/make_manifest.jl`
を書いたときに実際に踏んだ)。Windows のパスを値に入れれば `\\U` で同じく壊れる。
出荷 JSON の現行フィールドはたまたま無傷だったが、**公開する codec としては欠陥**。

RFC 8259 の必須分だけを扱う: `"` `\\` と制御文字 (U+0000..U+001F)。
非 ASCII はそのまま通す (UTF-8 のまま出すのが JSON として正しい)。

⚠ **出荷テーブルのバイトは変わらない。**エスケープが要る文字を含む文字列は
現行の一式に 1 つも無いことを実測済 (`tools/json_escape_audit.jl`)。"""
function json_escape(s::AbstractString)
    # 速い道: エスケープが要らなければ元の文字列をそのまま返す (割り当てゼロ)
    needs = false
    for c in s
        if c == '"' || c == '\\' || c < ' '
            needs = true
            break
        end
    end
    needs || return s
    out = IOBuffer()
    for c in s
        if c == '"'
            print(out, "\\\"")
        elseif c == '\\'
            print(out, "\\\\")
        elseif c == '\n'
            print(out, "\\n")
        elseif c == '\r'
            print(out, "\\r")
        elseif c == '\t'
            print(out, "\\t")
        elseif c == '\b'
            print(out, "\\b")
        elseif c == '\f'
            print(out, "\\f")
        elseif c < ' '
            print(out, "\\u", lpad(string(UInt16(c); base=16), 4, '0'))
        else
            print(out, c)
        end
    end
    return String(take!(out))
end

function write_json(io::IO, v; indent=0)
    pad = "  "^indent
    if v isa Dict
        println(io, "{")
        ks = sort(collect(keys(v)))
        for (i, k) in enumerate(ks)
            print(io, pad, "  \"", json_escape(string(k)), "\": ")
            write_json(io, v[k]; indent=indent + 1)
            println(io, i < length(ks) ? "," : "")
        end
        print(io, pad, "}")
    elseif v isa AbstractMatrix
        # 260806Cl 追加 (GOS 面): 行の配列として書く。行 = 第 1 添字なので
        # JSON 側は v[i][j] = v[i, j]。改行を入れるのは行の境目だけ
        println(io, "[")
        for i in axes(v, 1)
            print(io, pad, "  ")
            write_json(io, view(v, i, :); indent=indent + 1)
            println(io, i < last(axes(v, 1)) ? "," : "")
        end
        print(io, pad, "]")
    elseif v isa AbstractVector
        print(io, "[")
        for (i, x) in enumerate(v)
            write_json(io, x; indent=indent)
            i < length(v) && print(io, ", ")
        end
        print(io, "]")
    elseif v isa AbstractString
        print(io, "\"", json_escape(v), "\"")
    elseif v isa Bool || v isa Integer
        print(io, v)
    elseif v === nothing
        print(io, "null")
    else
        print(io, repr(Float64(v)))
    end
end
