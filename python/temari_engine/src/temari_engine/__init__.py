"""temari_engine — Temari のエンジンを Python から呼び、出力と出荷物を所属と役割を確かめて読む (L-C の E2、作者決定 I65)

- `run(args)` : Julia の CLI を subprocess で走らせ、envelope つきの出力を厳密に読む
- `read_output(path)` : 単発の出力 (`--json`) を読む (envelope 仕様 v1 の検査。研究入口の読み方)
- `load(path)` : 通常入口。一式の manifest で所属と `artifact_role` を確かめてから読む
  (`control`・所属の無いファイル・名乗りだけの `computed` は拒否。`research=True` で研究入口)
- `xray_lines.line_response(z, e0_keV, dataset_dir)` : dataset F v7.0.0 の殻ごとの σ と F(s) に xraylib の空孔の連鎖・蛍光収率・
  放射遷移の割合を掛けた、線ごとの X 線の応答 (作者決定 I83。xraylib は実行時の任意の依存)

依存は標準ライブラリだけ (xray_lines だけ任意で xraylib)。仕様 = Temari の `docs/notes/envelope_spec_v1_2026-09-22.md`。
"""
from ._strictjson import StrictJSONError, load_strict, loads_strict
from .envelope import Output, read_output, validate_envelope
from .errors import (EngineRunError, EnvelopeError, LineResponseError, MembershipError, RoleError, TemariEngineError,
                     XraylibError)
from .run import RunResult, default_repo, run
from .sets import (ArtifactSet, Member, files_digest, known_sets, load, load_f_channel, load_factors_release, load_set,
                   manifest_digest, media_types_digest, verify_f_release)

__version__ = "0.3.0"   # 0.3.0 (2026-09-25、作者決定 I83・I84): 線ごとの X 線の応答 xray_lines (xraylib は任意の依存)。
#   0.2.0 (2026-09-23、作者決定 I75): temari.artifact_set v2 を読む。
#   同日 I76 で v2 の digest に media_type を足した (0.2.0 は未公開なので版は据え置き)

__all__ = [
    "ArtifactSet", "EngineRunError", "EnvelopeError", "LineResponseError", "Member", "MembershipError", "Output", "RoleError",
    "RunResult", "XraylibError",
    "StrictJSONError", "TemariEngineError", "default_repo", "files_digest", "known_sets", "load", "load_f_channel",
    "load_factors_release", "load_set", "load_strict", "loads_strict", "manifest_digest", "media_types_digest", "read_output",
    "run",
    "validate_envelope",
    "verify_f_release",
]
