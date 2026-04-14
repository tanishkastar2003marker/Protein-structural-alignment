"""
Final robust implementation of τ-equivalence algorithm
Fully working + saves all results

"""

import numpy as np
from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.spatial.distance import cdist
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')


# ─────────────────────────────────────────────
# PART 1: Load protein structures
# ─────────────────────────────────────────────

def get_ca_atoms(pdb_file, chain_id=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_file)

    coords = []
    labels = []

    selected_chain = None
    first_chain = None

    for model in structure:
        for chain in model:
            print(f"DEBUG: {pdb_file} → found chain: {repr(chain.id)}")

            if first_chain is None:
                first_chain = chain

            if chain_id is not None and chain.id == chain_id:
                selected_chain = chain
                break
        break

    # ✅ fallback if chain not found
    if selected_chain is None:
        print(f"⚠️ Chain '{chain_id}' not found in {pdb_file}, using first chain")
        selected_chain = first_chain

    for res in selected_chain:
        if 'CA' in res:
            coords.append(res['CA'].get_vector().get_array())
            labels.append(f"{res.resname}{res.id[1]}")

    coords = np.array(coords)

    if coords.size == 0:
        print(f"❌ ERROR: No CA atoms in {pdb_file}")

    return coords, labels


# ─────────────────────────────────────────────
# PART 2: Distance + neighbors
# ─────────────────────────────────────────────

def build_distance_matrix(coords):
    if coords.ndim != 2 or len(coords) == 0:
        raise ValueError("Invalid coordinate array")
    return cdist(coords, coords)


def get_neighbors(dist_matrix, lo=5.0, hi=8.0):
    neighbors = {}
    n = len(dist_matrix)

    for i in range(n):
        neighbors[i] = [
            j for j in range(n)
            if i != j and lo <= dist_matrix[i][j] <= hi
        ]
    return neighbors


# ─────────────────────────────────────────────
# PART 3: τ-equivalence
# ─────────────────────────────────────────────

def is_tau_compatible(i, j, dist_a, dist_b, equiv_set, tau):
    for (eq_a, eq_b) in equiv_set:
        if abs(dist_a[i][eq_a] - dist_b[j][eq_b]) >= tau:
            return False
    return True


def expand_from_seed(seed_i, seed_j, dist_a, dist_b,
                     neighbors_a, neighbors_b, tau):

    equiv = [(seed_i, seed_j)]
    queue = [(seed_i, seed_j)]
    used_a = {seed_i}
    used_b = {seed_j}

    while queue:
        cur_a, cur_b = queue.pop(0)

        for nb_a in neighbors_a[cur_a]:
            if nb_a in used_a:
                continue

            for nb_b in neighbors_b[cur_b]:
                if nb_b in used_b:
                    continue

                if is_tau_compatible(nb_a, nb_b,
                                     dist_a, dist_b,
                                     equiv, tau):
                    equiv.append((nb_a, nb_b))
                    queue.append((nb_a, nb_b))
                    used_a.add(nb_a)
                    used_b.add(nb_b)
                    break

    return equiv


# ─────────────────────────────────────────────
# PART 4: Find equivalences
# ─────────────────────────────────────────────

def rank_by_degree(neighbors, ntop=10):
    return sorted(neighbors.keys(),
                  key=lambda x: len(neighbors[x]),
                  reverse=True)[:ntop]


def find_tau_equivalences(coords_a, coords_b, tau=1.5, ntop=10):

    dist_a = build_distance_matrix(coords_a)
    dist_b = build_distance_matrix(coords_b)

    neighbors_a = get_neighbors(dist_a)
    neighbors_b = get_neighbors(dist_b)

    seeds_a = rank_by_degree(neighbors_a, ntop)
    seeds_b = rank_by_degree(neighbors_b, ntop)

    best = []

    for i in seeds_a:
        for j in seeds_b:
            result = expand_from_seed(
                i, j, dist_a, dist_b,
                neighbors_a, neighbors_b, tau
            )
            if len(result) > len(best):
                best = result

    return best


# ─────────────────────────────────────────────
# PART 5: δ-equivalences
# ─────────────────────────────────────────────

def find_delta_equivalences(coords_a, coords_b,
                           tau_equiv, delta=3.0):

    if len(tau_equiv) < 3:
        return [], 999.0

    idx_a = [e[0] for e in tau_equiv]
    idx_b = [e[1] for e in tau_equiv]

    sup = SVDSuperimposer()
    sup.set(coords_a[idx_a], coords_b[idx_b])
    sup.run()

    rms = sup.get_rms()

    rot, tran = sup.get_rotran()
    coords_b_transformed = np.dot(coords_b, rot) + tran

    delta_equiv = []
    for ia, ib in zip(idx_a, idx_b):
        diff = np.linalg.norm(coords_a[ia] - coords_b_transformed[ib])
        if diff < delta:
            delta_equiv.append((ia, ib, round(diff, 2)))

    return delta_equiv, round(rms, 2)


# ─────────────────────────────────────────────
# PART 6: ANALYSIS (FINAL FIXED)
# ─────────────────────────────────────────────

def analyze_pair(pdb_a, pdb_b, name_a, name_b,
                 chain_a=None, chain_b=None,
                 tau_values=[0.5,1.0,1.5,2.0,2.5,3.0]):

    print(f"\n{'='*55}")
    print(f"{name_a}  vs  {name_b}")
    print(f"{'='*55}")

    coords_a, _ = get_ca_atoms(pdb_a, chain_a)
    coords_b, _ = get_ca_atoms(pdb_b, chain_b)

    print(f"Residues: {name_a}={len(coords_a)}, {name_b}={len(coords_b)}")

    if len(coords_a) == 0 or len(coords_b) == 0:
        print("❌ Skipping due to missing data")
        return None

    rows = []

    for tau in tau_values:
        print(f"Running τ = {tau} Å ...", end=' ')

        tau_eq = find_tau_equivalences(coords_a, coords_b, tau)
        delta_eq, rms = find_delta_equivalences(coords_a, coords_b, tau_eq)

        print(f"τ={len(tau_eq)}, δ={len(delta_eq)}, rms={rms}")

        rows.append({
            'τ (Å)': tau,
            'τ-equivalences': len(tau_eq),
            'δ-equivalences': len(delta_eq),
            'δ-r.m.s. (Å)': rms
        })

    df = pd.DataFrame(rows)

    # ✅ SAVE FILE (CRITICAL FIX)
    os.makedirs("results", exist_ok=True)
    safe_name = f"{name_a}_vs_{name_b}".replace(" ", "_")
    file_path = f"results/{safe_name}.csv"
    df.to_csv(file_path, index=False)

    print("\nResults Table:")
    print(df.to_string(index=False))
    print(f"✅ Saved: {file_path}")

    return df


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":

    pairs = [
        ("pdb_files/3dfr.pdb", "pdb_files/4dfr.pdb",
         "Lcasei DHFR", "Ecoli DHFR", "A", "A"),

        ("pdb_files/2hhb.pdb", "pdb_files/1mbd.pdb",
         "Haemoglobin-a", "Myoglobin", "A", "A"),

        ("pdb_files/2lzm.pdb", "pdb_files/1lyz.pdb",
         "T4 Lysozyme", "Hen Lysozyme", "A", "A"),

        ("pdb_files/3cyt.pdb", "pdb_files/3c2c.pdb",
         "Cytochrome-c", "Cytochrome-c2", None, None),

        ("pdb_files/2aza.pdb", "pdb_files/1pcy.pdb",
         "Azurin", "Plastocyanin", None, None),
    ]

    all_results = []

    for pdb_a, pdb_b, na, nb, ca, cb in pairs:
        df = analyze_pair(pdb_a, pdb_b, na, nb, ca, cb)

        if df is not None:
            row = df[df['τ (Å)'] == 1.5].iloc[0]
            all_results.append({
                'Protein Pair': f"{na} vs {nb}",
                'τ-equiv': row['τ-equivalences'],
                'δ-equiv': row['δ-equivalences'],
                'RMSD': row['δ-r.m.s. (Å)']
            })

    summary_df = pd.DataFrame(all_results)

    print("\n" + "="*60)
    print("FINAL SUMMARY (τ = 1.5 Å)")
    print("="*60)
    print(summary_df.to_string(index=False))

