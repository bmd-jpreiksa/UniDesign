from __future__ import annotations

import os
from pathlib import Path

from unidesign import DesignDomain, DesignProtein, Structure
from unidesign.api.design import SiteSpec


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    os.environ.setdefault("UNIDESIGN_LIBRARY_PATH", str(repo_root / "library"))

    pdb_path = repo_root / "example" / "MonomerDesign" / "1agy" / "1agy.pdb"

    structure = Structure.from_pdb(pdb_path)

    # Restrict a handful of design sites to specific amino-acid choices while
    # allowing other sites to sample the full alphabet.
    domain = DesignDomain(
        design_sites=[
            SiteSpec(chain="A", position=35, allowed="AVIL"),  # hydrophobic palette
            SiteSpec(chain="A", position=38, allowed="FWY"),   # aromatic-only
            SiteSpec(chain="A", position=71),                  # unrestricted
            SiteSpec(chain="A", position=90),                  # unrestricted
        ],
        repack_sites=[
            SiteSpec(chain="A", position=72),
            SiteSpec(chain="A", position=73),
        ],
    )

    designer = DesignProtein(structure, domain)
    designer.run(trajectories=1)

    print("Total energy:", designer.best_sequence_energy["total"])
    print("Sequence identity:", designer.best_sequence_energy["sequence_identity"])
    print("Restricted site choices:")
    for site in [("A", 35), ("A", 38), ("A", 71), ("A", 90)]:
        aa = designer.best_sequence[site]
        print(f"  {site[0]}:{site[1]:4d} -> {aa}")


if __name__ == "__main__":
    main()

