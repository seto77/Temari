# Temari 次フェーズ指示書 (2026-08-07)

*P3 (出口の拡張) 完了 + 原子場の物理の作り直し途中での引き継ぎ。
正本は `計画書.md` と `CLAUDE.md`、本書はそこからの差分・残務・作業順。
前回の引き継ぎは `docs/next_phase_2026-08-06.md` (P1 完了時点)。次セッションの冒頭で読むこと。*

---

## 1. 現在地

| 項目 | 状態 |
|---|---|
| リポジトリ | `https://github.com/seto77/Temari` (**近く public 化予定**)。CI + GitHub Pages 稼働中 |
| 出口 | **5 つ**: F(s,E0) / dσ/dΔE (EELS) / δ_l (弾性位相) / GOS / f_x・f_e (原子散乱因子) |
| 原子場 | **完全 Dirac SCF (DHFS) が既定**。`--nodscf` で旧非相対論 SCF |
| 交換 | **Xα (α=1) のまま**。厳密交換 (KLI) の機構は完成したが**まだ SCF に配線していない** |
| 検証 | selftest T0–T16 ALL PASS (~28 s) / refcheck WORST 9.044e-08 / カーネル === / 5 チャネル === |

### 今セッションで入ったもの (新しい順)

| コミット | 内容 |
|---|---|
| `e722f3b` | **整数占有の自己項補正**。これで全元素の交換漸近が −1/r に |
| `33c689b` | KLI 補正項 + 「平均配置の壁」の発見 |
| `3f9071e` | 平均配置の Slater ポテンシャル (角度係数を閉殻漸近で固定) |
| `b2e44c7`, `8bbb148` | 動径 Slater 関数 Y^k と自己相互作用の厳密相殺 |
| `9dae051` | **交換の診断書** (`docs/exchange_diagnosis_2026-08-07.md`) |
| `e64cc38` | 交換係数を引数化 + 自分で作ったバグ 2 件と既存バグ 1 件の修正 |
| `b6698a8`, `7a3de21`, `75537df` | Dirac SCF の実装 → 既定化 |
| `a63af77`, `5eb1608`, `4535ab5`, `fa1eeeb` | 出口 4 つ (f_x / GOS / δ_l / EELS) |
| `382f11a`, `a7499fe` | L0–L5 層分割 |
| `7501f01`, `e992376` | 球ベッセル Miller 規格化の 0/0 修正 (handout へも `0933c0e`) |

---

## 2. ★次の一手 — KLI を SCF へ配線し、Latter を削除する

**これだけやれば、この作業の目的が達成される。**機構は全部できていて、残るのは配線。

### やること

`SCFAtom` に `exchange::Symbol = :xalpha | :kli` を足し、`:kli` のとき:

```julia
vh = hartree(r, rho)
vx = kli_exchange_potential(P, q, l, r, eps)[1]   # 前反復の軌道から
veff = @. -z / r + vh + vx                        # ★min(...) の Latter クリップを外す
pot  = RvSpline(r, veff .* r, -(z - nel + 1))     # 漸近は物理が出す
```

**漸近が自動で正しくなる**のが要点: V_H → N/r、V_x → −1/r なので
V → −(Z−N+1)/r。中性なら −1/r、core-hole イオン (N=Z−1) なら −2/r で、
いま `latter_charge` に手で入れている値と厳密に一致する。

### 設計上の注意 (検討済み)

- **鶏と卵**: KLI は軌道を要るので、初回反復は Xα でブートストラップし、
  2 反復目から KLI に切り替える
- **混合**: 密度だけでなく**ポテンシャル vx も混合**する
  (`vx = (1-beta)*vx + beta*vx_new`)。KLI は密度からは決まらないため
- **Dirac 経路は別途**。`orbs[key]` に入っているのは占有加重和 (診断用) で
  軌道そのものではない。まず**非相対論経路で配線・検証**し、Dirac は次段
  (小成分 F を交換に入れるかの判断も要る)
- **コスト**: 副殻ペア × 多重極ごとに `ykr` を呼ぶ。Au (非相対論 14 副殻) で
  1 反復あたり ~1000 回 × 20000 点。10–30 s/原子を見込む

### 検査 (この順に)

1. **Latter 無しで収束すること**、かつ**有効ポテンシャルの漸近が −1/r** (中性)
2. **閉殻元素で Latter 版と大きく変わらないこと** (Ne, Ar は元々漸近が正しかった)
3. **開殻元素で密度が変わること** (C, Si, Au。ここが効くはず)
4. **payoff の測定**: f_x を公開パラメータ化と比べ、**α という knob 無しで**
   α=0.75 相当 (C 2.36% / Fe 0.43%) に迫るか。σ_own/σ_Bote も同時に測る
5. 既存ゲート一式 (下記 §5)

---

## 3. 今セッションで確定した物理 (蒸し返さない)

### 3.1 交換の診断 — 正本は `docs/exchange_diagnosis_2026-08-07.md`

- f_x の残差は**常に正**、絶対値の山は s ≈ 0.2–0.4 → **密度が収縮しすぎ**
- **Latter は主犯ではない** (Fe で有無が 1.55% vs 1.56%)。主因は局所交換の強さ
- **α を動かすのは解決にならない**: 最適 α は全元素 0.70–0.75 (Schwarz の値域) だが、
  f_x は 0.75 を、σ と 1s 固有値は 1 を好む。α=1+Latter は「固有値を束縛エネルギーに
  合わせる」設計、α≈0.7 は「HF の密度を再現する」値で、目標が違う

### 3.2 厳密交換の定式 (実装済み・検証済み)

    E_x = −(1/2) Σ_{ab} Σ_k c^k(l_a,l_b) W^k_{ab} G^k(ab)
    c^k = [3j(l_a k l_b;000)]²,  G^k(ab) = ∫ P_a P_b Y^k(ab;r)/r dr
    W^k_{ab} = q_a q_b/2 + δ_{ab}[(2l_a+1)q_a − q_a²/2]/(2k+1)

第 2 項が**整数占有の自己項補正**。球平均のために m の**選び方**を平均すると
⟨n(m)n(m′)⟩ = f² + δ_{mm′}(f−f²) となり、⟨n²⟩=⟨n⟩ (占有は整数) だから出る項。
角度因子は D_k(l) = Σ_m[3j(l k l;−m,0,m)]² = **1/(2k+1)** (l に依らない、実測確認)。

**この 1 項で 3 つが同時に正しくなる**:
1 電子の完全相殺 (4.4e-16) / 閉殻は補正ゼロ / **開殻の漸近が −1/r**
(C 2p²: −0.335 → −1.002、Au 6s¹: → −1.00001)。**スピン分極は不要**。

⚠ `33c689b` のコミットメッセージに「球平均した開殻では −1/r は原理的に回復
できない」と書いたが**これは誤り**。`e722f3b` で訂正済み。

### 3.3 危険な落とし穴 (実際に踏んだ)

- **u_x の 1/2**: δE/δP_a = 2q_aε_aP_a なので u_x P_a = (1/(2q_a))δE_x/δP_a。
  落とすと Slater ポテンシャルと 2 倍食い違う。恒等式 **Σ_a q_a ū_a = 2E_x** で固定 (T16)
- **u_x 自体を作らない**: P_a の節で 0/0。常に **P_a² u_x** の形で扱う
- **交換係数の二重計上**: `slater_vx` に α を畳み込むと終状態の KS(2/3) が (2/3)α に
  なる。**`slater_vx` は素の Slater 形**、α は SCF 側、2/3 は終状態側で掛ける (T13b)
- **キャッシュキーの取り違え**: 束縛状態の鍵に SCF 種別・α を入れ忘れると、
  `--dscf` 実行が入れた 1s を非 dscf 実行が黙って読む。refcheck が
  `|dE_b/E_b| = 2.2e-3` で検出した (K 殻だけ、L 殻無傷という症状)

---

## 4. その先の残務 (優先度順)

### 4.1 Dirac 経路への KLI 拡張

非相対論で効果が出たら。小成分 F を交換に入れるか (Dirac–Fock 本来は入る) の
判断が要る。Au 1s では ∫F² が全体の ~9%。

### 4.2 スピン分極 SCF (診断書の案 B)

作者了解済み「必ずやることになる」。**ただし漸近の問題は §3.2 で解けたので、
当初想定より優先度は下がった。**残る動機は開殻の交換エネルギー自体の精度
(Hund 則の効果) と、将来の磁性・スピン分解量。

### 4.3 P3 の残り / P4

- **d²σ/dΩdΔE** (小)。K=0 分岐が既に θ グリッド上で S/Q⁴ を評価しているので配管のみ
- **σ(β, Δ)** (中)。EELS 定量の k-factor。ロードマップ最大の裾野
- **Mott 弾性断面積** (中)。`phase` 出口の δ_l から。ただし**スピン分解が要る**
  (κ 分解 Dirac 連続状態) のと、低 l の主値を連続化する必要がある (Levinson)
- **イオンの f_x** (小〜中)。`SCFAtom` は任意占有を受けるので陽イオンは半日。
  陰イオンは安定化モデル (Watson 球等) を**明示的に宣言**してから

### 4.4 v4 物理 + M 殻

正本は ReciPro 側 `.project-guidance/ReciPro/ReciPro_STEM-EDX_v4精度検討.md`。
−Re(DX*) 干渉項が残っている最大の物理項。

### 4.5 E6 (R テーブル共有) / 未検証の高速化 16 案

`docs/speedup_audit_2026-08-05.json` の verdict 無し項目。v4 の生成時間に直結。

---

## 5. 作業の掟 (前回から継続 + 今回の追加)

- **Julia は `julia +1.11` (1.11.9) が最終ゲート**。PATH の `julia` は 1.12.6
- **検証手順**:
  ```powershell
  julia +1.11 -t 4 src/ionization.jl selftest      # T0-T16、~28 s
  julia +1.11 -t 4 src/ionization.jl refcheck      # WORST 9.044e-08 が基準値
  julia +1.11 -t 1 tools/verify_simd_bessel.jl     # 288 ケース
  julia +1.11 -t 1 tools/verify_e5_qlane.jl        #  75 ケース
  julia +1.11 -t 4 tools/bitident_snapshot.jl before.txt   # ★変更の「前」に取る
  ```
- **refcheck は `dirac_scf=false, x_alpha=1.0` に固定してある**。実装間の比較であり、
  Python 参照値は非相対論・α=1 の処方で作られているため。**処方を変えてもここは動かさない**
- **値が変わる修正では「変化を無効化した版」も走らせる** (例: `J0_MIN=0.0`)。
  それが修正前とビット同一なら、差分は全て意図した変化に帰属できる
- **処方を変えたら model_id に印を付ける** — `model_id_of(rel, dscf, x_alpha)` が唯一の
  組み立て口。表示も JSON も必ずここを通す (分岐ごとに文字列連結を書き足して事故った)
- **キャッシュキーには処方の全次元を入れる** (SCF 種別・α・スキーマ版)。
  入れ忘れると黙って別処方の原子を読む
- **公開パラメータ化の係数はリポに入れない**。比較スクリプトは scratchpad のみ。
  参照値は `C:\Users\seto\source\repos\Crystallography\Crystallography\Atom\` を使う
  (作者指示。Cromer–Mann / Waasmaier–Kirfel / Peng / Kirkland)
- コミットメッセージは英語 + `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`

---

## 6. コマンド早見表

```powershell
# 出口
julia +1.11 -t 4 src/ionization.jl      26 K 200 --rel        # F(s) (EDX)
julia +1.11 -t 4 src/ionization.jl edge 26 K 200 --rel        # dσ/dΔE (EELS)
julia +1.11 -t 4 src/ionization.jl gos  26 K --epsmax 2000    # GOS (E0 を取らない)
julia +1.11 -t 4 src/ionization.jl phase 26 100               # δ_l (Z と ε[eV])
julia +1.11 -t 4 src/ionization.jl fx   26                    # f_x / f_e (Z のみ)

# 処方の切り替え
--nodscf     非相対論 SCF に戻す (既定は完全 Dirac SCF)
--rel        放出電子のスカラー相対論 (v3 処方)
--nonrel     fx 出口のみ: 非相対論密度で比較
```

---

## 7. 参照

- **交換の診断書**: `docs/exchange_diagnosis_2026-08-07.md` (次の一手の根拠)
- 前回の引き継ぎ: `docs/next_phase_2026-08-06.md`
- 正本: `計画書.md`、`docs/architecture.md` (層とファイルの対応)
- 検証の一覧: `docs/src/en/verification.md` (T0–T16 の期待値つき)
- 来歴と運用注意: `src/IMPORT.md`、`CLAUDE.md`
