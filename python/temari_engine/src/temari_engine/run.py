"""Julia のエンジン (`src/ionization.jl`) を subprocess で 1 回走らせ、envelope つきの出力を厳密に読む (E2)"""
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass

from .envelope import Output, read_output
from .errors import EngineRunError, EnvelopeError


def default_repo():
    """環境変数 TEMARI_REPO、無ければこのファイルから上へ辿って `src/ionization.jl` のある dir"""
    env = os.environ.get("TEMARI_REPO")
    if env:
        return os.path.abspath(env)
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.isfile(os.path.join(d, "src", "ionization.jl")):
            return d
        d = os.path.dirname(d)
    raise EngineRunError("Temari の repo が見つからない (TEMARI_REPO を指定する)")


@dataclass(frozen=True)
class RunResult:
    output: Output
    returncode: int
    stdout: str
    stderr: str
    argv: tuple


def run(args, repo=None, julia="julia", threads=4, workdir=None, timeout=None, allow_nonzero=False):
    """`julia src/ionization.jl <args…> --json <一時ファイル>` を走らせて RunResult を返す。

    - `args` は出口の引数 (例: `["mott", "79", "10000"]`、F(s) は `["26", "K", "200", "--quick"]`)。`--json` は渡さない。
    - `julia` は実行ファイル名か、`["julia", "+1.11.9"]` のような列。
    - `workdir` は Julia の cwd (atom_cache は cwd に対する相対パス)。既定は repo の根 (CLI の文書と同じ)。
    - 出力の envelope の `command` が `args + ["--json", "<json>"]` と一致することを確かめる (別の木・別の入口を走らせていない)。
    - 非 0 終了は既定で EngineRunError (mott の `truncated` は exit 2 で JSON を書く — `allow_nonzero=True` なら結果を返す)。
    """
    args = [str(a) for a in args]
    if "--json" in args:
        raise ValueError("args に --json を入れない (run が一時ファイルへ書かせる)")
    repo = os.path.abspath(repo) if repo else default_repo()
    entry = os.path.join(repo, "src", "ionization.jl")
    if not os.path.isfile(entry):
        raise EngineRunError("エンジンの入口が無い: %s" % entry)
    jl = [julia] if isinstance(julia, str) else list(julia)
    tmp = tempfile.mkdtemp(prefix="temari_engine_")
    out = os.path.join(tmp, "out.json")
    argv = tuple(jl + ["--startup-file=no", "--project=%s" % repo, "-t", str(threads), entry] + args + ["--json", out])
    try:
        try:
            p = subprocess.run(argv, cwd=workdir or repo, capture_output=True, text=True, encoding="utf-8",
                               errors="replace", timeout=timeout)
        except FileNotFoundError as e:
            raise EngineRunError("julia を起動できない: %s" % e)
        except subprocess.TimeoutExpired as e:
            raise EngineRunError("時間切れ (%s 秒)" % timeout, stdout=e.stdout or "", stderr=e.stderr or "")
        if not os.path.isfile(out):
            raise EngineRunError("出力が書かれなかった (exit %d)" % p.returncode, p.returncode, p.stdout, p.stderr)
        try:
            o = read_output(out)
        except EnvelopeError as e:
            raise EngineRunError("出力の envelope が不正 (exit %d): %s" % (p.returncode, e), p.returncode, p.stdout, p.stderr)
        want = args + ["--json", "<json>"]
        if o.envelope["command"] != want:
            raise EngineRunError("出力の command %r が呼び出し %r と違う" % (o.envelope["command"], want),
                                 p.returncode, p.stdout, p.stderr)
        o = Output(payload=o.payload, envelope=o.envelope, path=None)   # 一時ファイルは消す
        if p.returncode != 0 and not allow_nonzero:
            raise EngineRunError("エンジンが非 0 で終了 (exit %d)" % p.returncode, p.returncode, p.stdout, p.stderr)
        return RunResult(output=o, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr, argv=argv)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
