#!/usr/bin/env python3
"""Prepare a deterministic 5-target/100-ligand DUD-E benchmark."""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Lipinski
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmark-100"
SOURCE = BENCH / "source" / "diverse"
PREPARED = BENCH / "prepared"
TARGETS = ("ampc", "cxcr4", "gcr", "hivpr", "hivrt")
LIGANDS_PER_TARGET = 20


def infer_element(line: str) -> str:
    atom_name = line[12:16].strip()
    residue = line[17:20].strip()
    letters = "".join(c for c in atom_name if c.isalpha())
    if line.startswith("ATOM"):
        return (letters[:1] or "C").upper()
    two_letter = {"CL", "BR", "ZN", "MG", "MN", "FE", "CA", "NA", "CU", "CO", "NI", "SE"}
    candidate = letters[:2].upper()
    if candidate in two_letter and (residue == candidate or atom_name.upper().startswith(candidate)):
        return candidate.title()
    return (letters[:1] or "C").upper()


def normalize_pdb(source: Path, destination: Path) -> None:
    residue_aliases = {"ALB": "ALA", "ASM": "ASN", "ASQ": "ASP", "GLZ": "GLY", "SEM": "SER", "TYM": "TYR"}
    lines = []
    for raw in source.read_text().splitlines():
        if not raw.startswith(("ATOM", "HETATM")):
            lines.append(raw)
            continue
        line = raw.ljust(80)
        residue = line[17:20].strip()
        if residue == "WAU":
            continue
        if residue in residue_aliases:
            line = line[:17] + residue_aliases[residue].rjust(3) + line[20:]
        element = infer_element(line)
        lines.append(line[:76] + element.rjust(2) + line[78:])
    destination.write_text("\n".join(lines) + "\n")


def prepare_receptor(target: str, target_dir: Path) -> tuple[Path, Path]:
    normalized = target_dir / "receptor_normalized.pdb"
    normalize_pdb(SOURCE / target / "receptor.pdb", normalized)
    receptor = target_dir / "receptor.pdbqt"
    command = ["obabel", str(normalized), "-O", str(receptor), "-xr"]
    completed = subprocess.run(command, text=True, capture_output=True)
    (target_dir / "receptor_prepare.log").write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"receptor preparation failed for {target}; see {target_dir/'receptor_prepare.log'}")
    crystal = Chem.MolFromMol2File(str(SOURCE / target / "crystal_ligand.mol2"), removeHs=False, sanitize=False)
    if crystal is None:
        raise RuntimeError(f"could not parse crystal ligand for {target}")
    positions = crystal.GetConformer().GetPositions()
    minimum, maximum = positions.min(axis=0), positions.max(axis=0)
    center, size = (minimum + maximum) / 2, (maximum - minimum) + 12.0
    box = target_dir / "receptor.box.txt"
    box.write_text("\n".join([
        f"center_x = {center[0]:.3f}", f"center_y = {center[1]:.3f}", f"center_z = {center[2]:.3f}",
        f"size_x = {size[0]:.3f}", f"size_y = {size[1]:.3f}", f"size_z = {size[2]:.3f}",
    ]) + "\n")
    return receptor, box


def prepare_ligands(target: str, target_dir: Path) -> list[dict]:
    ligand_dir = target_dir / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.ForwardSDMolSupplier(gzip.open(SOURCE / target / "actives_final.sdf.gz", "rb"), removeHs=False)
    preparator = MoleculePreparation()
    selected = []
    for source_index, molecule in enumerate(supplier):
        if molecule is None:
            continue
        heavy_atoms = molecule.GetNumHeavyAtoms()
        rotors = Lipinski.NumRotatableBonds(molecule)
        if heavy_atoms > 64 or rotors > 8:
            continue
        name = molecule.GetProp("_Name").strip() or f"active_{source_index:04d}"
        safe_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        try:
            setups = preparator.prepare(molecule)
            if len(setups) != 1:
                continue
            pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
            if not ok:
                continue
        except Exception:
            continue
        branches = sum(line.startswith("BRANCH") for line in pdbqt.splitlines())
        atom_types = sorted({line.split()[-1] for line in pdbqt.splitlines() if line.startswith(("ATOM", "HETATM"))})
        unsupported = any(atom_type.startswith("G") or atom_type == "W" for atom_type in atom_types)
        if branches > 8 or unsupported:
            continue
        output = ligand_dir / f"{len(selected)+1:02d}_{safe_name}.pdbqt"
        output.write_text(pdbqt)
        selected.append({
            "target": target, "ligand": safe_name, "source_index": source_index,
            "pdbqt": str(output.relative_to(ROOT)), "heavy_atoms": heavy_atoms,
            "rotatable_bonds": rotors, "branches": branches, "atom_types": atom_types,
        })
        if len(selected) == LIGANDS_PER_TARGET:
            break
    if len(selected) != LIGANDS_PER_TARGET:
        raise RuntimeError(f"only found {len(selected)} supported ligands for {target}")
    return selected


def write_maps(target: str, receptor: Path, box: Path, ligands: list[dict], target_dir: Path) -> Path:
    maps_dir = target_dir / "maps"
    maps_dir.mkdir(exist_ok=True)
    prefix = maps_dir / "affinity"
    logs = []
    for entry in ligands:
        command = [str(ROOT / "bin/vina"), "--receptor", str(receptor),
                   "--ligand", str(ROOT / entry["pdbqt"]), "--config", str(box),
                   "--write_maps", str(prefix), "--force_even_voxels", "--score_only", "--verbosity", "0"]
        completed = subprocess.run(command, text=True, capture_output=True)
        logs.append(completed.stdout + completed.stderr)
        # DUD-E's stored conformer may be outside the crystal-ligand box. Vina writes
        # the maps before score_only checks that pose, so validate the actual artifact.
        if not any(maps_dir.glob("affinity.*.map")):
            (target_dir / "map_generation.log").write_text("\n".join(logs))
            raise RuntimeError(f"map generation failed for {target}")
    (target_dir / "map_generation.log").write_text("\n".join(logs))
    return maps_dir


def main() -> None:
    PREPARED.mkdir(parents=True, exist_ok=True)
    manifest = {"benchmark_version": "dude-diverse-100-v1", "source": "DUD-E Diverse",
                "targets": list(TARGETS), "ligands_per_target": LIGANDS_PER_TARGET, "cases": []}
    for target in TARGETS:
        print(f"Preparing {target}", flush=True)
        target_dir = PREPARED / target
        target_dir.mkdir(parents=True, exist_ok=True)
        receptor, box = prepare_receptor(target, target_dir)
        ligands = prepare_ligands(target, target_dir)
        maps = write_maps(target, receptor, box, ligands, target_dir)
        for entry in ligands:
            entry.update({"receptor": str(receptor.relative_to(ROOT)), "config": str(box.relative_to(ROOT)),
                          "maps": str(maps.relative_to(ROOT))})
        manifest["cases"].extend(ligands)
    (BENCH / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {len(manifest['cases'])} docking cases across {len(TARGETS)} proteins")


if __name__ == "__main__":
    main()
