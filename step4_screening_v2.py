"""
Step 4 スクリーニングパイプライン v2
修正内容：
  - コンセンサル残基スコアを熱安定性軸から削除
  - コンセンサルを「5/5完全保存位置のハードフィルタ」に移動
  - 熱安定性軸に表面電荷スコア（TGP論文知見）を追加

使い方:
    python step4_screening_v2.py \
        --fasta mpnn_output/seqs/2B3P.fa \
        --pdb 2B3P.pdb \
        --thermompnn_dir ./ThermoMPNN \
        --output_dir ./screening_results
"""

import os, sys, csv, json, argparse, subprocess
import numpy as np
from collections import Counter
from Bio import Align

# ─── 参照配列 ─────────────────────────────────────────────────────────
GFP_SEQS = {
    "sfGFP":   "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "avGFP":   "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "amacGFP": "MSKGEELFTGIVPVLIELDGDVHGHKFSVRGEGEGDADYGKLEIKFICTTGKLPVPWPTLVTTLSYGILCFARYPEHMKMNDFFKSAMPEGYIQERTIFFQDDGKYKTRGEVKFEGDTLVNRIELKGMDFKEDGNILGHKLEYNFNSHNVYIMPDKANNGLKVNFKIRHNIEGGGVQLADHYQTNVPLGDGPVLIPINHYLSCQTAISKDRNETRDHMVFLEFFSACGHTHGMDELYK",
    "cgreGFP": "MTALTEGAKLFEKEIPYITELEGDVEGMKFIIKGEGTGDATTGTIKAKYICTTGDLPVPWATILSSLSYGVFCFAKYPRHIADFFKSTQPDGYSQDRIISFDNDGQYDVKAKVTYENGTLYNRVTVKGTGFKSNGNILGMRVLYHSPPHAVYILPDRKNGGMKIEYNKAFDVMGGGHQMARHAQFNKPLGAWEEDYPLYHHLTVWTSFGKDPDDDETDHLTIVEVIKAVDLETYR",
    "ppluGFP": "MPAMKIECRITGTLNGVEFELVGGGEGTPEQGRMTNKMKSTKGALTFSPYLLSHVMGYGFYHFGTYPSGYENPFLHAINNGGYTNTRIEKYEDGGVLHVSFSYRYEAGRVIGDFKVVGTGFPEDSVIFTDKIIRSNATVEHLHPMGDNVLVGSFARTFSLRDGGYYSFVVDSHMHFKSAIHPSILQNGGPMFAFRRVEELHSNTELGIVEYQHAFKTPIAFA",
}
SFGFP = GFP_SEQS["sfGFP"]

# sfGFP の固定すべき残基（1-indexed）- ProteinMPNN側で固定済み
FIXED_POSITIONS = {
    65: "T", 66: "Y", 67: "G",  # クロモフォアトライアド
    30: "R",                     # S30R イオンペアネットワーク
    99: "S",                     # F99S cycle-3
    153: "T",                    # M153T cycle-3
    163: "A",                    # V163A cycle-3
}

MIN_LEN, MAX_LEN = 220, 250

# sfGFP の net charge ベースライン
SFGFP_NET_CHARGE = (SFGFP.count('D') + SFGFP.count('E')) \
                 - (SFGFP.count('K') + SFGFP.count('R'))  # = +6（正電荷過多）


# ═══════════════════════════════════════════════════════════════
# 起動時に一度だけ計算：5/5完全保存位置テーブル
# ═══════════════════════════════════════════════════════════════

def build_perfect_conserved_positions() -> dict:
    """
    5種GFPのアライメントから「5/5完全保存位置」を抽出
    → ハードフィルタとして使用（熱安定性スコアとしては使わない）
    
    注意：クロモフォア65番はsfGFP固有のS65T変異のため5/5保存ではない
          → 既存のFIXED_POSITIONSで別途管理
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    position_aas = {i: [] for i in range(len(SFGFP))}
    for name, seq in GFP_SEQS.items():
        best = aligner.align(SFGFP, seq)[0]
        sfgfp_pos = 0
        for s_aa, t_aa in zip(*best):
            if s_aa != "-" and t_aa != "-":
                position_aas[sfgfp_pos].append(t_aa)
                sfgfp_pos += 1
            elif s_aa != "-":
                sfgfp_pos += 1

    perfect = {
        pos: aas[0]
        for pos, aas in position_aas.items()
        if len(aas) == 5 and len(set(aas)) == 1
    }
    return perfect  # {0-indexed pos: aa}


print("Building conservation table...")
PERFECT_CONSERVED = build_perfect_conserved_positions()
print(f"  5/5 perfectly conserved positions: {len(PERFECT_CONSERVED)}")


# ═══════════════════════════════════════════════════════════════
# A. ハードフィルタ（v2：コンセンサル保存位置チェック追加）
# ═══════════════════════════════════════════════════════════════

def parse_fasta(fasta_path: str) -> list:
    seqs = []
    header, seq = None, ""
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq:
                    seqs.append({"header": header, "seq": seq})
                header, seq = line[1:], ""
            else:
                seq += line
    if header and seq:
        seqs.append({"header": header, "seq": seq})
    print(f"\n[Input] {len(seqs)} sequences loaded")
    return seqs


def extract_mpnn_score(header: str) -> float:
    for tok in header.split(","):
        tok = tok.strip()
        if tok.startswith("score="):
            try:
                return float(tok.split("=")[1])
            except ValueError:
                pass
    return 999.0


def hard_filter(seqs: list) -> list:
    """
    ハードフィルタ（v2）
    既存5条件 + 新規：5/5完全保存位置チェック

    5/5完全保存位置でコンセンサルと異なる残基を持つ配列は排除
    （機能・折り畳みに必須の残基が変異している → 輝度ゼロリスク）
    """
    passed = []
    reasons = Counter()

    for entry in seqs:
        seq = entry["seq"].replace("/", "").replace("-", "")

        # 1. 長さ
        if not (MIN_LEN <= len(seq) <= MAX_LEN):
            reasons["length"] += 1; continue

        # 2. Mスタート
        if not seq.startswith("M"):
            reasons["no_M"] += 1; continue

        # 3. クロモフォアトライアド TYG
        if seq[64:67] != "TYG":
            reasons["chromophore"] += 1; continue

        # 4. sfGFP 安定化変異4点
        fail = False
        for pos1, aa in FIXED_POSITIONS.items():
            if len(seq) >= pos1 and seq[pos1-1] != aa:
                reasons[f"stab_mutation"] += 1
                fail = True; break
        if fail:
            continue

        # 5. WTと完全一致除外
        if seq == SFGFP:
            reasons["identical_WT"] += 1; continue

        # 6. ★新規：5/5完全保存位置チェック
        #    クロモフォア・安定化変異は上記で確認済みなのでスキップ
        conserved_fail = False
        for pos0, cons_aa in PERFECT_CONSERVED.items():
            if pos0 < len(seq) and seq[pos0] != cons_aa:
                # 既存チェック済み位置は除外
                if (pos0+1) not in FIXED_POSITIONS:
                    reasons["conserved_violation"] += 1
                    conserved_fail = True
                    break
        if conserved_fail:
            continue

        entry["seq"] = seq
        entry["mpnn_score"] = extract_mpnn_score(entry["header"])
        entry["identity"] = sum(a==b for a,b in zip(seq, SFGFP)) / len(SFGFP)
        passed.append(entry)

    print(f"\n[A] Hard filter: {len(seqs)} → {len(passed)}")
    print(f"    Rejection reasons: {dict(reasons)}")
    return passed


# ═══════════════════════════════════════════════════════════════
# B. ESM-1v スコア（輝度制約）
# ═══════════════════════════════════════════════════════════════

def compute_esm1v_scores(seqs: list, batch_size: int = 16) -> list:
    try:
        import torch
        import esm as esm_lib
    except ImportError:
        print("[B] ESM not available. Skipping.")
        for e in seqs:
            e["esm1v_score"] = 0.0
            e["esm1v_vs_wt"] = 0.0
        return seqs

    print("\n[B] Loading ESM-1v...")
    model, alphabet = esm_lib.pretrained.esm1v_t33_650M_UR90S_1()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)

    def score_batch(batch_seqs):
        data = [(f"s{i}", s) for i, s in enumerate(batch_seqs)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            logits = model(tokens, repr_layers=[], return_contacts=False)["logits"]
            log_probs = torch.log_softmax(logits, dim=-1)
        scores = []
        for b, seq in enumerate(batch_seqs):
            ll = sum(log_probs[b, i+1, alphabet.get_idx(aa)].item()
                     for i, aa in enumerate(seq))
            scores.append(ll / len(seq))
        return scores

    wt_score = score_batch([SFGFP])[0]
    print(f"    sfGFP WT baseline: {wt_score:.4f}")

    all_seqs = [e["seq"] for e in seqs]
    all_scores = []
    for i in range(0, len(all_seqs), batch_size):
        all_scores.extend(score_batch(all_seqs[i:i+batch_size]))

    for entry, score in zip(seqs, all_scores):
        entry["esm1v_score"] = score
        entry["esm1v_vs_wt"] = score - wt_score

    threshold = wt_score - 0.5
    filtered = [e for e in seqs if e["esm1v_score"] >= threshold]
    print(f"\n[B] ESM-1v filter: {len(seqs)} → {len(filtered)}")
    return filtered


def filter_by_mpnn_score(seqs: list, top_n: int = 50) -> list:
    ranked = sorted(seqs, key=lambda x: x.get("mpnn_score", 999))[:top_n]
    print(f"\n[B'] MPNN score filter: {len(seqs)} → {len(ranked)}")
    return ranked


# ═══════════════════════════════════════════════════════════════
# C. 熱安定性スコアリング（v2）
#    ThermoMPNN ΔΔG + 表面電荷スコア
#    ※コンセンサルスコアはハードフィルタに移動済み
# ═══════════════════════════════════════════════════════════════

def run_thermompnn_ssm(pdb_path: str, thermompnn_dir: str, out_dir: str) -> dict:
    """ThermoMPNN SSM を sfGFP 構造に実行 → ΔΔG テーブルを取得"""
    ssm_csv = os.path.join(out_dir, "sfgfp_ssm.csv")

    if not os.path.exists(thermompnn_dir):
        print(f"\n[C] ThermoMPNN not found. Using dummy ΔΔG=0.")
        return {}

    if not os.path.exists(ssm_csv):
        print(f"\n[C] Running ThermoMPNN SSM on sfGFP...")
        script = os.path.join(thermompnn_dir, "analysis", "custom_inference.py")
        cmd = [sys.executable, script,
               "--pdb", pdb_path, "--chain", "A", "--out", ssm_csv]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    Failed: {result.stderr[:200]}")
            return {}
    else:
        print(f"\n[C] Using cached ThermoMPNN SSM: {ssm_csv}")

    ssm_table = {}
    try:
        with open(ssm_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                pos = int(row.get("position", row.get("pos", 0)))
                mut_aa = row.get("MUT_AA", row.get("mut", ""))
                ddg = float(row.get("ddG", row.get("ddg", 0.0)))
                ssm_table[(pos, mut_aa)] = ddg
        print(f"    Loaded {len(ssm_table)} entries")
    except Exception as e:
        print(f"    Parse error: {e}")
    return ssm_table


def compute_thermompnn_ddg(seq: str, ssm_table: dict) -> float:
    """
    生成配列の累積 ΔΔG（加法近似）
    sfGFP と異なる位置のみ ΔΔG を加算
    負値 = 安定化（良い）
    """
    if not ssm_table:
        return 0.0
    return sum(
        ssm_table.get((i+1, s_aa), 0.0)
        for i, (s_aa, ref_aa) in enumerate(zip(seq, SFGFP))
        if s_aa != ref_aa
    )


def compute_charge_score(seq: str) -> float:
    """
    表面電荷スコア（TGP論文知見ベース）

    TGP の設計思想：
      - 正電荷（K/R）を減らし負電荷（D/E）を増やす
      - eCGP123 の pI≈7.0 → TGP の net charge = -10
      - 凝集抑制・熱安定性向上に寄与

    sfGFP の net charge = +6（正電荷過多）
    生成配列が sfGFP より負に動いているほどスコア高

    スコア = (seq_net_charge - sfgfp_net_charge) の符号を反転
            = sfgfp_net_charge - seq_net_charge
            → 正値 = sfGFP より負電荷が増えている（良い）
    """
    seq_charge = (seq.count('D') + seq.count('E')) \
               - (seq.count('K') + seq.count('R'))
    # sfGFP より負になっているほど正のスコア
    return float(SFGFP_NET_CHARGE - seq_charge)


def compute_thermal_scores(seqs: list, ssm_table: dict) -> list:
    """全配列に ThermoMPNN ΔΔG と表面電荷スコアを計算"""
    for e in seqs:
        e["thermompnn_ddg"]  = compute_thermompnn_ddg(e["seq"], ssm_table)
        e["charge_score"]    = compute_charge_score(e["seq"])
        e["net_charge"]      = (e["seq"].count('D') + e["seq"].count('E')) \
                             - (e["seq"].count('K') + e["seq"].count('R'))

    # sfGFP ベースラインを表示
    sfgfp_charge = compute_charge_score(SFGFP)
    print(f"\n[C] Thermal scoring complete")
    print(f"    sfGFP net charge: {SFGFP_NET_CHARGE} (baseline=0)")
    print(f"    charge_score > 0 = more negative than sfGFP (good)")
    return seqs


# ═══════════════════════════════════════════════════════════════
# D. 複合スコアランキング（v2）
# ═══════════════════════════════════════════════════════════════

def rank_candidates(seqs: list, top_n: int = 15) -> list:
    """
    複合スコア（v2）
    輝度制約：ESM-1v
    熱安定性：ThermoMPNN ΔΔG × 表面電荷スコア

    composite = (esm1v_normalized)
              × exp(-ddg * 0.1)          # ddg負=安定化ボーナス
              × sigmoid(charge_score)    # 電荷改善ボーナス
    """
    import math

    def sigmoid(x):
        return 1 / (1 + math.exp(-x * 0.3))

    for e in seqs:
        esm_norm = max(0.0, e.get("esm1v_vs_wt", 0.0) + 1.0)
        ddg_bonus = math.exp(-e.get("thermompnn_ddg", 0.0) * 0.1)
        charge_bonus = sigmoid(e.get("charge_score", 0.0))
        e["composite_score"] = esm_norm * ddg_bonus * charge_bonus

    ranked = sorted(seqs, key=lambda x: -x["composite_score"])[:top_n]

    print(f"\n[D] Final ranking: top {len(ranked)} candidates")
    print(f"\n{'Rank':>4} {'MPNN':>7} {'ESM-1v':>8} "
          f"{'ΔΔG':>7} {'NetChg':>7} {'Comp':>7} {'ID':>5}")
    print("-" * 52)
    for i, e in enumerate(ranked):
        print(f"{i+1:>4} {e.get('mpnn_score',0):>7.4f} "
              f"{e.get('esm1v_score',0):>8.4f} "
              f"{e.get('thermompnn_ddg',0):>7.3f} "
              f"{e.get('net_charge',0):>7} "
              f"{e['composite_score']:>7.4f} "
              f"{e['identity']:>5.3f}")
    return ranked


# ═══════════════════════════════════════════════════════════════
# メイン
# ═══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",          required=True)
    parser.add_argument("--pdb",            required=True)
    parser.add_argument("--thermompnn_dir", default="./ThermoMPNN")
    parser.add_argument("--output_dir",     default="./screening_results_v2")
    parser.add_argument("--top_n",          type=int, default=15)
    parser.add_argument("--esm_top_n",      type=int, default=50)
    parser.add_argument("--skip_esm",       action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # A. ハードフィルタ（コンセンサル保存位置チェック込み）
    raw      = parse_fasta(args.fasta)
    filtered = hard_filter(raw)

    # B. ESM-1v + MPNN スコア（輝度制約）
    if not args.skip_esm:
        filtered = compute_esm1v_scores(filtered)
    filtered = filter_by_mpnn_score(filtered, top_n=args.esm_top_n)

    # C. 熱安定性スコアリング（ThermoMPNN + 表面電荷）
    ssm_table = run_thermompnn_ssm(args.pdb, args.thermompnn_dir, args.output_dir)
    filtered  = compute_thermal_scores(filtered, ssm_table)

    # D. 複合スコアランキング
    final = rank_candidates(filtered, top_n=args.top_n)

    # 保存
    out_json = os.path.join(args.output_dir, "screening_results_v2.json")
    with open(out_json, "w") as f:
        json.dump([{k:v for k,v in e.items() if k!="header"}
                   for e in final], f, indent=2)

    out_csv = os.path.join(args.output_dir, "candidates_for_md.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank","seq","mpnn_score","esm1v_score",
                    "thermompnn_ddg","net_charge","charge_score",
                    "composite_score","identity"])
        for i, e in enumerate(final):
            w.writerow([i+1, e["seq"],
                        e.get("mpnn_score",""), e.get("esm1v_score",""),
                        e.get("thermompnn_ddg",""), e.get("net_charge",""),
                        e.get("charge_score",""), e.get("composite_score",""),
                        e.get("identity","")])

    print(f"\n✓ Saved: {out_json}")
    print(f"✓ Saved: {out_csv}")
    print(f"\n→ 次ステップ: candidates_for_md.csv → ESMFold → MD @ 72℃")


if __name__ == "__main__":
    main()
