# L0 JSON — 標準ライブラリに JSON が無いための最小パーサ / ライタ
#
# 用途は bote_salvat.json・reference_values.json の読み込みと --json 保存のみ。
# 数値ではないが、同じく「Python の標準品を自前で置き換えた小道具」なので L0 に置く。

# ---- 最小 JSON パーサ (bote_salvat.json / reference_values.json 専用) ----
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
                i += 1
                ch = b[i]
                write(buf, ch == UInt8('n') ? 0x0a :
                           ch == UInt8('t') ? 0x09 : ch)
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

"JSON ファイルを読む (このパッケージの 2 つの JSON 専用の最小実装)"
parse_json_file(path::String) = _json_value(read(path), 1)[1]

# ---- 最小 JSON writer (--json 保存用) ----
function write_json(io::IO, v; indent=0)
    pad = "  "^indent
    if v isa Dict
        println(io, "{")
        ks = sort(collect(keys(v)))
        for (i, k) in enumerate(ks)
            print(io, pad, "  \"", k, "\": ")
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
        print(io, "\"", v, "\"")
    elseif v isa Bool || v isa Integer
        print(io, v)
    elseif v === nothing
        print(io, "null")
    else
        print(io, repr(Float64(v)))
    end
end
