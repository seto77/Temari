"""temari_engine の例外。どれも TemariEngineError の子"""


class TemariEngineError(Exception):
    """temari_engine が出す例外の親"""


class EnvelopeError(TemariEngineError):
    """単発の出力の envelope が仕様 v1 に合わない (欠落・知らない版・知らない欄・exit の食い違い など)"""


class MembershipError(TemariEngineError):
    """所属を確かめられない (manifest が無い・既知の一式でない・sha256 や digest が合わない・manifest に載っていない)"""


class RoleError(TemariEngineError):
    """通常入口が読まない役割 (control・既知の表に無い computed・所属の無い単発の出力)"""


class XraylibError(TemariEngineError):
    """xraylib (Python の束縛) が無い、または要るデータを返さない (線ごとの応答 `xray_lines` の任意の依存)"""


class LineResponseError(TemariEngineError):
    """線ごとの応答を出せない (E₀ が表の行に無い・対応しない線・strict で供給元の殻が無い など)"""


class EngineRunError(TemariEngineError):
    """Julia のエンジンの実行が失敗した (起動できない・非 0 終了・出力が無い / 壊れている)"""

    def __init__(self, message, returncode=None, stdout="", stderr=""):
        super().__init__(message)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class FetchError(TemariEngineError):
    """公開物を取得できない (どの写しからも取れない・大きさや sha256 が台帳と合わない・HTTPS でない URL など)。`temari_engine.data`"""


class OfflineError(TemariEngineError):
    """offline の指定 (引数 `offline=True` か環境変数 `TEMARI_OFFLINE=1`) で、要る公開物が cache に無い。`temari_engine.data`"""


class ArchiveError(TemariEngineError):
    """書庫を安全に展開できない (危ない要素・想定外の形・展開の失敗)。`temari_engine.data`"""


class SpeciesLookupError(TemariEngineError):
    """元素と電荷などから表の種を 1 つに決められない (候補が無い・2 つ以上ある)。`candidates` に候補の一覧を持つ"""

    def __init__(self, message, candidates=()):
        super().__init__(message)
        self.candidates = tuple(candidates)
