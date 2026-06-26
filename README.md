# GFP Design Pipeline — 2026 SynBio Challenges

Computational pipeline for designing novel GFP variants with high initial brightness and extreme thermal stability (72 °C resistance), submitted to the 2026 Protein Design in SynBio Challenges.

---

## Strategy

The competition score simplifies to:

$$\text{Score} = \frac{F_\text{initial}}{F_\text{initial,WT}} \times \frac{F_\text{final}}{F_\text{initial}} = \frac{F_\text{final}}{F_\text{initial,WT}}$$

Since $F_\text{initial}$ cancels, **thermal stability is the dominant driver**. We split the 6 submission slots into two portfolios:

| Slots | Approach | Goal |
|-------|----------|------|
| 1–3 | Conservative: sfGFP + ProteinMPNN | High-probability stable & bright |
| 4–6 | Exploratory: TGP-style surface charge redesign | Maximum thermal stability |

---

## Environment Setup

```bash
# Create isolated conda environment
conda create -n gfp_design python=3.10 -y
conda activate gfp_design

# PyTorch (CUDA 12.x)
pip install torch --index-url https://download.pytorch.org/whl/cu128

# Core dependencies
pip install numpy scipy biopython fair-esm \
            pytorch-lightning omegaconf tqdm \
            ml-collections pandas wandb transformers

# OpenMM (CUDA-enabled, via conda-forge)
conda install -c conda-forge openmm -y

# Clone external tools (place alongside gfp_design/)
git clone https://github.com/dauparas/ProteinMPNN.git   ../ProteinMPNN
git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git ../ThermoMPNN

# Patch ThermoMPNN to use local ProteinMPNN weights
sed -i 's|model_weight_dir = os.path.join(cfg.platform.thermompnn_dir.*|model_weight_dir = "../ProteinMPNN/vanilla_model_weights"|' \
    ../ThermoMPNN/transfer_model.py
```

---

## Pipeline Overview

```
sfGFP structure (PDB: 2B3P)
        │
        ▼
[Phase 1] ProteinMPNN — inverse folding
        │   6,000 candidate sequences
        ▼
[Phase 2] Screening — 6,000 → 15 candidates
        │   A. Hard filter (length / M-start / chromophore / stabilizing mutations)
        │   B. Brightness constraint (ESM-1v log-likelihood + ProteinMPNN score)
        │   C. Thermal stability scoring (ThermoMPNN SSM + surface charge)
        ▼
[Phase 3] MD @ 72 °C — 15 → 3 final sequences
            10 ns NPT simulation per candidate
            CA-RMSD vs. sfGFP WT baseline (mean = 8.40 Å)
```

---

## Running the Pipeline

### Step 1 — Download structure & generate sequences

```bash
cd gfp_design

# Download sfGFP crystal structure
curl "https://files.rcsb.org/download/2B3P.pdb" -o 00_setup/2B3P.pdb

# Run ProteinMPNN (≈ 2–3 min on RTX 4090 / Blackwell)
bash run_pipeline.sh
```

`run_pipeline.sh` fixes residues 30, 99, 153, 163, 65–67 (chromophore triad + sfGFP stabilizing mutations) and samples at temperatures 0.10 / 0.15 / 0.20, generating 2,000 sequences each.

### Step 2 — Post-process & screen

```bash
# Fix missing Met1 and restore chromophore TYG
python3 fix_fasta.py \
    --input  01_proteinmpnn/seqs/2B3P.fa \
    --output 01_proteinmpnn/seqs/2B3P_fixed.fa

# Run ThermoMPNN SSM on sfGFP (once, ~5 s)
python ../ThermoMPNN/analysis/custom_inference.py \
    --pdb        00_setup/2B3P.pdb \
    --chain      A \
    --model_path ../ThermoMPNN/models/thermoMPNN_default.pt \
    --out_dir    02_screening/

# Screen candidates
python step4_screening_v2.py \
    --fasta          01_proteinmpnn/seqs/2B3P_fixed.fa \
    --pdb            00_setup/2B3P.pdb \
    --thermompnn_dir ../ThermoMPNN \
    --output_dir     02_screening
```

### Step 3 — MD simulation @ 72 °C

```bash
# Run MD for top 15 candidates (≈ 30 min per sequence on Blackwell 6000)
python 03_md_simulation/run_md.py
```

Sequences with mean CA-RMSD below sfGFP WT baseline (8.40 Å) are selected.

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Start from sfGFP (2B3P) | Competition WT baseline; already highly stable (ΔG = 9.57 kcal/mol) |
| Fix chromophore TYG | Fluorescence requires intact T65-Y66-G67 |
| Fix S30R | 5-residue ion-pair network; primary sfGFP stabilizer (Pédelacq et al. 2006) |
| Fix F99S / M153T / V163A | cycle-3 stabilizing mutations |
| ThermoMPNN as filter, not sum | Model is designed for single point mutations; summing 70+ ΔΔGs violates additivity |
| MD @ 72 °C as final filter | Directly measures kinetic stability; matches competition assay condition |

---

## MD Protocol

| Parameter | Value |
|-----------|-------|
| Force field | AMBER14 + TIP3P-FB |
| System size | ~49,000 atoms |
| Temperature | 345.15 K (72 °C) |
| Ensemble | NPT |
| Time step | 2 fs |
| Production | 10 ns |
| Metric | CA-RMSD vs. initial structure |

sfGFP WT baseline: **mean CA-RMSD = 8.40 Å** (same protocol, same GPU).

---

## Final Results (Slots 1–3)

| Slot | Sequence ID | mean CA-RMSD | vs. WT |
|------|-------------|-------------|--------|
| 1 | rank9 | 6.79 Å | −19% |
| 2 | seq2  | 6.87 Å | −18% |
| 3 | rank3 | 7.36 Å | −12% |

Submitted sequences: `results/submission.csv`

---

## Repository Structure

```
gfp_design/
├── 00_setup/                  sfGFP structure (2B3P.pdb)
├── 01_proteinmpnn/            ProteinMPNN outputs (6,000 sequences)
├── 02_screening/              Screening results & ThermoMPNN SSM table
├── 03_md_simulation/          MD scripts & RMSD results
├── 04_slots4_6_tgp/           Exploratory design (TGP-style, in progress)
├── 05_final_submission/       submission.csv
├── results/                   submission.csv (competition format)
├── run_pipeline.sh            ProteinMPNN execution
├── step4_screening_v2.py      Screening pipeline
└── README.md
```

---

## References

1. Pédelacq et al. *Nature Biotechnology* (2006) — sfGFP, stabilizing mutations
2. Sarkisyan et al. *Nature* (2016) — avGFP fitness landscape (training data)
3. Close et al. *Proteins* (2015) — TGP surface charge engineering
4. Dauparas et al. — ProteinMPNN
5. Kuhlman Lab — ThermoMPNN
