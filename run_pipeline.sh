#!/bin/bash
set -e

PYTHON="/home/student/miniconda3/envs/gfp_design/bin/python"
WORKDIR="$HOME/gfp_design"
MPNN_DIR="$HOME/ProteinMPNN"
PDB_DIR="$WORKDIR/00_setup"
OUT_DIR="$WORKDIR/01_proteinmpnn"
FIXED="$OUT_DIR/fixed_positions.jsonl"
PARSED="$OUT_DIR/parsed_pdbs.jsonl"

echo "=== Step 1: PDB 確認 ==="
ls -lh "$PDB_DIR/2B3P.pdb"

echo ""
echo "=== Step 2: PDB → JSONL 変換 ==="
$PYTHON "$MPNN_DIR/helper_scripts/parse_multiple_chains.py" \
    --input_path "$PDB_DIR" \
    --output_path "$PARSED"

echo ""
echo "=== Step 3: 固定位置 JSONL 作成 ==="
$PYTHON -c "
import json, sys
fixed = {'2B3P': {'A': [65, 66, 67, 30, 99, 153, 163]}}
with open('$FIXED', 'w') as f:
    f.write(json.dumps(fixed) + '\n')
print('Written:', '$FIXED')
"

echo ""
echo "=== Step 4: ProteinMPNN 実行 ==="
$PYTHON "$MPNN_DIR/protein_mpnn_run.py" \
    --jsonl_path            "$PARSED" \
    --fixed_positions_jsonl "$FIXED" \
    --out_folder            "$OUT_DIR" \
    --num_seq_per_target    2000 \
    --sampling_temp         "0.10 0.15 0.20" \
    --use_soluble_model \
    --seed                  42 \
    --batch_size            16 \
    --omit_AAs              "X"

echo ""
echo "=== 完了 ==="
echo "出力: $OUT_DIR/seqs/2B3P.fa"
