"""
MD @ 72℃ - 2B3P起点版
候補配列の変異を 2B3P 結晶構造に適用してから MD 実行
"""

import os, csv, sys
import numpy as np
from multiprocessing import Process


def apply_mutations_and_run_md(ref_pdb, seq, seq_name, gpu_id, out_dir, n_steps=5_000_000):
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer

    os.makedirs(out_dir, exist_ok=True)

    SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

    AA3 = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS",
           "I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN","P":"PRO","Q":"GLN",
           "R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR"}

    print(f"[{seq_name}] GPU {gpu_id} 開始", flush=True)

    # ── 1. 変異リスト作成（2B3P は残基2スタート）──
    mutations = []
    for i in range(1, min(len(seq), len(SFGFP))):
        if seq[i] != SFGFP[i]:
            pdb_resnum = i + 1
            mutations.append(f"{AA3[SFGFP[i]]}-{pdb_resnum}-{AA3[seq[i]]}")

    print(f"[{seq_name}] 変異数: {len(mutations)}", flush=True)

    # ── 2. pdbfixer で構造準備 ──
    # CRO/ACY等の非標準残基を除去してから変異適用
    fixer = PDBFixer(filename=ref_pdb)
    fixer.removeHeterogens(False)   # CRO, ACY 等を全除去
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    fixer.applyMutations(mutations, "A")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    fixer.addSolvent(
        padding=1.2 * unit.nanometer,
        ionicStrength=0.15 * unit.molar
    )
    print(f"[{seq_name}] 構造準備完了", flush=True)

    # ── 3. 力場・システム ──
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
    )

    T = 345.15 * unit.kelvin  # 72℃
    integrator = mm.LangevinMiddleIntegrator(T, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
    system.addForce(mm.MonteCarloBarostat(1.0*unit.atmospheres, T))

    platform = mm.Platform.getPlatformByName("CUDA")
    props = {"CudaDeviceIndex": str(gpu_id), "CudaPrecision": "mixed"}
    simulation = app.Simulation(fixer.topology, system, integrator, platform, props)
    simulation.context.setPositions(fixer.positions)

    # ── 4. エネルギー最小化 ──
    print(f"[{seq_name}] エネルギー最小化...", flush=True)
    simulation.minimizeEnergy(maxIterations=1000)

    # ── 5. 加熱（300 → 345 K）──
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    for temp in range(300, 346, 5):
        integrator.setTemperature(temp * unit.kelvin)
        simulation.step(2500)
    print(f"[{seq_name}] 加熱完了", flush=True)

    # ── 6. 平衡化（0.5 ns）──
    simulation.step(250_000)
    print(f"[{seq_name}] 平衡化完了", flush=True)

    # ── 7. 本番MD（10 ns）──
    state0 = simulation.context.getState(getPositions=True)
    pos0 = np.array(state0.getPositions().value_in_unit(unit.nanometer))
    ca_indices = [a.index for a in simulation.topology.atoms() if a.name == "CA"]

    simulation.reporters.append(app.DCDReporter(f"{out_dir}/{seq_name}.dcd", 50_000))
    simulation.reporters.append(app.StateDataReporter(
        f"{out_dir}/{seq_name}.log", 50_000,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True
    ))

    rmsd_data = []
    report_interval = 250_000  # 0.5 ns ごと
    print(f"[{seq_name}] 本番MD @ 72℃ 開始", flush=True)

    for i in range(n_steps // report_interval):
        simulation.step(report_interval)
        state = simulation.context.getState(getPositions=True)
        pos   = np.array(state.getPositions().value_in_unit(unit.nanometer))
        ca_pos0 = pos0[ca_indices]
        ca_pos  = pos[ca_indices]
        rmsd = np.sqrt(((ca_pos - ca_pos0)**2).sum(axis=1).mean()) * 10  # Å
        t_ns = (i+1) * report_interval * 2e-6
        rmsd_data.append({"time_ns": round(t_ns, 2), "rmsd_A": round(rmsd, 3)})
        print(f"[{seq_name}] t={t_ns:.1f} ns  CA-RMSD={rmsd:.2f} Å", flush=True)

    with open(f"{out_dir}/{seq_name}_rmsd.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["time_ns", "rmsd_A"])
        w.writeheader(); w.writerows(rmsd_data)

    mean_rmsd = np.mean([r["rmsd_A"] for r in rmsd_data])
    print(f"[{seq_name}] 完了: mean CA-RMSD={mean_rmsd:.2f} Å", flush=True)


if __name__ == "__main__":
    REF_PDB = "/home/student/gfp_design/00_setup/2B3P.pdb"
    OUT_BASE = "/home/student/gfp_design/03_md_simulation/md_results"

    seqs = {
        "seq1": "MGKGAELLEGEVPILVELEGDVNGHKFSIRGEGKGIAAEGKLELKFVCTTGKLPVPWPTLVTTLTYGINCFAKYPEHMQEHDFFKACLPEGYRRTLTLSFKDDGTFETEAEVRFEGDTLVNRIKLKGTGFKEGGNILGHKIKYTYESFTVNITADAAANGIKATFTLKLPLEDGSVQEVKVEGRYTPIGAGPATLPAPHYLKVERELSRDPNEKRDHMVLHEKITAGGIAAP",
        "seq2": "MGKGAELLKGVVPVRVELEGDVNGHKFSIRGEGEGDAEEGRLRLKFVCTTGKLPVPWPTLVTTLTYGLSCFAKYPEHMQDHDFFKACMPEGYRRERKLSFKDDGTYYTEAEVRFEGDTLVNRIKLEGVGFKEGGNILGHKLEYSYESFTVNITADAAANGITAKFTLKLPVKDGSTQLVDVEERNTPIGEGPATLPQPHYLKVEVKLSKDPNEKRDHMVLEEYVTAGGIAAP",
        "seq3": "MGKGDELLKGEVPLEVRLEGDFNGHKFSVRGEGKGDAEKGLQHLKFVCTTGKLPVPWPTLVTTLTYGVLCTAKYPEHMKDHDFFKACLPEGYIREQTLSFKDDGTYHVKARVYFEGDTLVNEIELKGTGFKEGGNILGHKLKYTYNSYTVNITADEKNNGIKATYTIELPVEDGSTQLVDHEGTYTPIGEGPDKLPEPHYLKVEVKLSKDPNEKRDHMVLERKVTAGGIAAP",
    }

    jobs = [
        (REF_PDB, seqs["seq1"], "seq1", 0, f"{OUT_BASE}/seq1"),
        (REF_PDB, seqs["seq2"], "seq2", 1, f"{OUT_BASE}/seq2"),
        (REF_PDB, seqs["seq3"], "seq3", 0, f"{OUT_BASE}/seq3"),
    ]

    # seq1 + seq2 並列
    p1 = Process(target=apply_mutations_and_run_md, args=jobs[0])
    p2 = Process(target=apply_mutations_and_run_md, args=jobs[1])
    p1.start(); p2.start()
    p1.join();  p2.join()

    # seq3
    p3 = Process(target=apply_mutations_and_run_md, args=jobs[2])
    p3.start(); p3.join()

    print("\n=== 全MD完了 ===")
