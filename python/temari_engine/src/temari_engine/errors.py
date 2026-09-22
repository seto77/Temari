"""temari_engine の例外。どれも TemariEngineError の子"""


class TemariEngineError(Exception):
    """temari_engine が出す例外の親"""


class EnvelopeError(TemariEngineError):
    """単発の出力の envelope が仕様 v1 に合わない (欠落・知らない版・知らない欄・exit の食い違い など)"""


class MembershipError(TemariEngineError):
    """所属を確かめられない (manifest が無い・既知の一式でない・sha256 や digest が合わない・manifest に載っていない)"""


class RoleError(TemariEngineError):
    """通常入口が読まない役割 (control・既知の表に無い computed・所属の無い単発の出力)"""


class EngineRunError(TemariEngineError):
    """Julia のエンジンの実行が失敗した (起動できない・非 0 終了・出力が無い / 壊れている)"""

    def __init__(self, message, returncode=None, stdout="", stderr=""):
        super().__init__(message)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
