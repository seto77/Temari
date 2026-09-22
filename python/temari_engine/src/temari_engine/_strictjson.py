"""厳密な JSON の読み手 (標準ライブラリだけ)。

拒否するもの: UTF-8 でないバイト・BOM・重複キー・`NaN` / `Infinity` の字句・binary64 に収まらない数 (`1e400`)・構文の誤り。
Temari の出力は非有限の値を `null` にして envelope の `nonfinite` に位置を書くので、厳密に読めないものは壊れている。
"""
import json
import math


class StrictJSONError(ValueError):
    """規格外の JSON。`kind` は encoding / bom / duplicate_key / constant / nonfinite / syntax"""

    def __init__(self, kind, message):
        super().__init__(message)
        self.kind = kind


def _pairs(pairs):
    d = {}
    for k, v in pairs:
        if k in d:
            raise StrictJSONError("duplicate_key", "重複キー %r" % k)
        d[k] = v
    return d


def _constant(name):
    raise StrictJSONError("constant", "JSON の規格外の定数 %s" % name)


def _float(token):
    x = float(token)
    if not math.isfinite(x):
        raise StrictJSONError("nonfinite", "binary64 に収まらない数 %s" % token)
    return x


def loads_strict(data):
    """bytes (または str) を厳密に読む"""
    if isinstance(data, (bytes, bytearray)):
        if bytes(data[:3]) == b"\xef\xbb\xbf":
            raise StrictJSONError("bom", "UTF-8 の BOM がある")
        try:
            text = bytes(data).decode("utf-8")
        except UnicodeDecodeError as e:
            raise StrictJSONError("encoding", "UTF-8 でない: %s" % e)
    else:
        text = data
    try:
        return json.loads(text, object_pairs_hook=_pairs, parse_constant=_constant, parse_float=_float)
    except StrictJSONError:
        raise
    except ValueError as e:
        raise StrictJSONError("syntax", "JSON の構文の誤り: %s" % e)


def load_strict(path):
    with open(path, "rb") as f:
        return loads_strict(f.read())
