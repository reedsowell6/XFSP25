# polymer_cli.py

import random as rnd
import statistics
from polymerClases import macroMolecule, Position

# Seed the global RNG once for reproducible results
rnd.seed(42)

def ask_int(prompt, default):
    """Ask for an integer, use default if blank or invalid."""
    s = input(f"{prompt} ({default})?: ").strip()
    if not s:
        return default
    try:
        return int(s)
    except ValueError:
        print(f"  → invalid, using {default}")
        return default

def main():
    # 1) get user inputs
    target_N = ask_int("degree of polymerization", 1000)
    num_molecules = ask_int("How many molecules", 50)

    # 2) run the simulations
    centers = []
    e2e = []
    rg = []
    Mi = []  # molecular weights for each chain

    for _ in range(num_molecules):
        chain = macroMolecule(target_N)
        chain.freelyJointedChainModel()
        centers.append(chain.centerOfMass)
        e2e.append(chain.endToEndDistance)
        rg.append(chain.radiusOfGyration)
        Mi.append(chain.N * chain.merWt)

    # 3) average center-of-mass in nm
    avg_com = Position(
        x=sum(p.x for p in centers) / num_molecules,
        y=sum(p.y for p in centers) / num_molecules,
        z=sum(p.z for p in centers) / num_molecules
    )
    avg_com_nm = Position(
        x=avg_com.x * 1e9,
        y=avg_com.y * 1e9,
        z=avg_com.z * 1e9
    )

    # 4) convert lengths to micrometers
    e2e_um = [d * 1e6 for d in e2e]
    rg_um  = [r * 1e6 for r in rg]

    # 5) compute statistics
    avg_e2e = statistics.mean(e2e_um)
    std_e2e = statistics.stdev(e2e_um)
    avg_rg  = statistics.mean(rg_um)
    std_rg  = statistics.stdev(rg_um)

    # 6) compute PDI = Mw/Mn
    Mn = statistics.mean(Mi)
    Mw = sum(m*m for m in Mi) / sum(Mi)
    pdi = Mw / Mn

    # 7) print results
    print(f"Metrics for {num_molecules} molecules of degree of polymerization = {target_N}")
    print(f"Avg. Center of Mass (nm) = {avg_com_nm.x:.3f}, {avg_com_nm.y:.3f}, {avg_com_nm.z:.3f}")
    print("End-to-end distance (μm):")
    print(f"    Average = {avg_e2e:.3f}")
    print(f"    Std. Dev. = {std_e2e:.3f}")
    print("Radius of gyration (μm):")
    print(f"    Average = {avg_rg:.3f}")
    print(f"    Std. Dev. = {std_rg:.3f}")
    print(f"PDI = {pdi:.2f}")

if __name__ == "__main__":
    main()
