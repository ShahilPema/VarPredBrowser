#!/usr/bin/env python3
"""
Fetch AlphaFold structure for a gene and extract pLDDT scores.

AlphaFold stores pLDDT (predicted Local Distance Difference Test) scores
in the B-factor column of the PDB file. Values range from 0-100:
- >90: Very high confidence
- 70-90: Confident
- 50-70: Low confidence
- <50: Very low confidence

Usage:
    python scripts/fetch_alphafold.py [--gene GENE] [--uniprot UNIPROT]
"""

import argparse
import json
import requests
from pathlib import Path
from typing import Dict

# Try to import from browser config, fall back to local path
try:
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from browser.backend.config import get_data_dir, get_structures_dir
    DATA_DIR = get_data_dir()
    STRUCTURES_DIR = get_structures_dir()
except ImportError:
    DATA_DIR = Path(__file__).parent.parent / 'data'
    STRUCTURES_DIR = DATA_DIR / 'structures'

# Default configuration (SCN2A)
DEFAULT_UNIPROT_ID = 'Q99250'
DEFAULT_GENE_SYMBOL = 'SCN2A'


def download_alphafold_structure(uniprot_id: str, output_dir: Path, version: int = 6) -> Path:
    """Download AlphaFold PDB file for a UniProt accession."""

    alphafold_id = f"AF-{uniprot_id}-F1"
    url = f"https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_v{version}.pdb"

    output_path = output_dir / f"{alphafold_id}-model_v{version}.pdb"

    if output_path.exists():
        print(f"Structure already exists: {output_path}")
        return output_path

    print(f"Downloading AlphaFold structure: {alphafold_id}")
    print(f"  URL: {url}")

    response = requests.get(url)
    response.raise_for_status()

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path.write_text(response.text)

    print(f"  Saved to: {output_path}")
    return output_path


def extract_plddt_scores(pdb_path: Path) -> Dict[int, float]:
    """
    Extract pLDDT scores from AlphaFold PDB file.

    In AlphaFold PDB files, pLDDT is stored in the B-factor column.
    We extract one score per residue (using CA atoms).
    """

    plddt_scores = {}

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()

                # Use CA (alpha carbon) for residue-level score
                if atom_name == 'CA':
                    residue_num = int(line[22:26].strip())
                    b_factor = float(line[60:66].strip())
                    plddt_scores[residue_num] = b_factor

    return plddt_scores


def fetch_alphafold_metadata(uniprot_id: str) -> Dict:
    """Fetch metadata from AlphaFold API."""

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

    print(f"Fetching AlphaFold metadata...")
    response = requests.get(url)
    response.raise_for_status()

    data = response.json()
    if isinstance(data, list) and len(data) > 0:
        return data[0]
    return data


def main(gene_symbol: str, uniprot_id: str):
    print("=" * 70)
    print(f"Fetching AlphaFold Structure for {gene_symbol} ({uniprot_id})")
    print("=" * 70)

    # Create directories
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    # Download structure
    print("\n[Step 1] Downloading AlphaFold structure...")
    pdb_path = download_alphafold_structure(uniprot_id, STRUCTURES_DIR)

    # Extract pLDDT scores
    print("\n[Step 2] Extracting pLDDT scores...")
    plddt_scores = extract_plddt_scores(pdb_path)

    print(f"  Residues with pLDDT: {len(plddt_scores)}")

    # Calculate statistics
    scores = list(plddt_scores.values())
    mean_plddt = sum(scores) / len(scores)
    min_plddt = min(scores)
    max_plddt = max(scores)

    # Count by confidence category
    very_high = sum(1 for s in scores if s > 90)
    confident = sum(1 for s in scores if 70 <= s <= 90)
    low = sum(1 for s in scores if 50 <= s < 70)
    very_low = sum(1 for s in scores if s < 50)

    print(f"\n  pLDDT Statistics:")
    print(f"    Mean: {mean_plddt:.1f}")
    print(f"    Min:  {min_plddt:.1f}")
    print(f"    Max:  {max_plddt:.1f}")
    print(f"\n  Confidence Distribution:")
    print(f"    Very high (>90):  {very_high:4d} ({100*very_high/len(scores):.1f}%)")
    print(f"    Confident (70-90): {confident:4d} ({100*confident/len(scores):.1f}%)")
    print(f"    Low (50-70):       {low:4d} ({100*low/len(scores):.1f}%)")
    print(f"    Very low (<50):    {very_low:4d} ({100*very_low/len(scores):.1f}%)")

    # Fetch metadata
    print("\n[Step 3] Fetching AlphaFold metadata...")
    try:
        metadata = fetch_alphafold_metadata(uniprot_id)
        print(f"  Model version: {metadata.get('latestVersion', 'N/A')}")
        print(f"  UniProt: {metadata.get('uniprotAccession', 'N/A')}")
    except Exception as e:
        print(f"  Warning: Could not fetch metadata: {e}")
        metadata = {}

    # Save pLDDT scores as JSON for easy access
    print("\n[Step 4] Saving pLDDT scores...")
    plddt_path = DATA_DIR / f'{gene_symbol.lower()}_plddt.json'
    with open(plddt_path, 'w') as f:
        json.dump({
            'gene_symbol': gene_symbol,
            'uniprot_accession': uniprot_id,
            'alphafold_id': f"AF-{uniprot_id}-F1",
            'protein_length': len(plddt_scores),
            'mean_plddt': round(mean_plddt, 2),
            'min_plddt': round(min_plddt, 2),
            'max_plddt': round(max_plddt, 2),
            'plddt_by_residue': plddt_scores
        }, f, indent=2)
    print(f"  Saved to: {plddt_path}")

    # Also update/create structure_metadata.json
    print("\n[Step 5] Updating structure metadata...")
    structure_metadata_path = DATA_DIR / 'structure_metadata.json'

    if structure_metadata_path.exists():
        with open(structure_metadata_path, 'r') as f:
            structure_metadata = json.load(f)
    else:
        structure_metadata = {}

    structure_metadata.update({
        'gene_symbol': gene_symbol,
        'uniprot_accession': uniprot_id,
        'protein_length': len(plddt_scores),
        'plddt_by_residue': plddt_scores,
        'structures': {
            'alphafold': {
                'id': f"AF-{uniprot_id}-F1",
                'version': metadata.get('latestVersion', 6),
                'file': f"structures/AF-{uniprot_id}-F1-model_v{metadata.get('latestVersion', 6)}.pdb",
                'mean_plddt': round(mean_plddt, 2),
                'min_plddt': round(min_plddt, 2),
                'max_plddt': round(max_plddt, 2),
                'coverage': 1.0,
                'confidence_distribution': {
                    'very_high': very_high,
                    'confident': confident,
                    'low': low,
                    'very_low': very_low
                }
            }
        }
    })

    with open(structure_metadata_path, 'w') as f:
        json.dump(structure_metadata, f, indent=2)
    print(f"  Saved to: {structure_metadata_path}")

    print("\n" + "=" * 70)
    print("AlphaFold structure fetch complete!")
    print("=" * 70)

    return pdb_path, plddt_scores


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch AlphaFold structure and extract pLDDT scores')
    parser.add_argument('--gene', default=DEFAULT_GENE_SYMBOL, help='Gene symbol')
    parser.add_argument('--uniprot', default=DEFAULT_UNIPROT_ID, help='UniProt accession')

    args = parser.parse_args()
    main(args.gene, args.uniprot)
