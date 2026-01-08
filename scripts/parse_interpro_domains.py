#!/usr/bin/env python3
"""
InterPro Domain Parser for MANE Select Transcripts

Parses match_complete.xml.gz from InterPro to extract domain annotations
for MANE Select proteins. Uses representative domains for Domain/Family/Repeat
types to avoid redundancy.

Data sources:
- InterPro: https://ftp.ebi.ac.uk/pub/databases/interpro/releases/
- MANE: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
- UniProt ID mapping: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

Usage:
    python parse_interpro_domains.py \
        --interpro-xml /path/to/match_complete.xml.gz \
        --mane-summary /path/to/MANE.GRCh38.summary.txt.gz \
        --uniprot-mapping /path/to/HUMAN_9606_idmapping.dat.gz \
        --output /path/to/mane_domains.parquet

    # Or with automatic downloads:
    python parse_interpro_domains.py \
        --data-dir /path/to/data \
        --output /path/to/mane_domains.parquet \
        --download

Based on parsing approach from: https://github.com/rcalef/magneton
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import requests
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('parse_interpro_domains.log')
    ]
)
logger = logging.getLogger(__name__)

# Download URLs
INTERPRO_XML_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/releases/103.0/match_complete.xml.gz"
MANE_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"
UNIPROT_MAPPING_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"

# InterPro types that use the representative field for de-duplication
REPRESENTATIVE_TYPES = {"Domain", "Family", "Repeat"}

# Types to exclude (Family is too broad/generic)
EXCLUDED_TYPES = {"Family"}


def download_file(url: str, output_path: Path, description: str = None) -> Path:
    """
    Download a file using aria2c (if available) or requests.

    Args:
        url: URL to download
        output_path: Path to save the file
        description: Description for progress bar

    Returns:
        Path to downloaded file
    """
    if output_path.exists():
        logger.info(f"File already exists: {output_path}")
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)
    description = description or output_path.name

    logger.info(f"Downloading {description}: {url}")

    # Try aria2c first (faster for large files)
    aria2c_path = shutil.which('aria2c')
    if aria2c_path:
        try:
            cmd = [
                'aria2c',
                '--max-connection-per-server=8',
                '--split=8',
                '--min-split-size=1M',
                '--file-allocation=none',
                '--continue=true',
                '--auto-file-renaming=false',
                '-d', str(output_path.parent),
                '-o', output_path.name,
                url
            ]
            logger.info(f"Using aria2c for faster download...")
            result = subprocess.run(cmd, capture_output=False)
            if result.returncode == 0 and output_path.exists():
                logger.info(f"Downloaded: {output_path}")
                return output_path
        except Exception as e:
            logger.warning(f"aria2c failed: {e}, falling back to requests")

    # Fall back to requests with progress bar
    response = requests.get(url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get('content-length', 0))

    with open(output_path, 'wb') as f:
        with tqdm(total=total_size, unit='iB', unit_scale=True, desc=description) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

    logger.info(f"Downloaded: {output_path}")
    return output_path


def ensure_data_files(data_dir: Path, download: bool = False) -> tuple[Path, Path, Path]:
    """
    Ensure all required data files exist, downloading if necessary.

    Args:
        data_dir: Directory to store data files
        download: Whether to download missing files

    Returns:
        Tuple of (interpro_xml_path, mane_summary_path, uniprot_mapping_path)
    """
    data_dir = Path(data_dir)
    raw_dir = data_dir / "raw"

    interpro_xml = raw_dir / "match_complete.xml.gz"
    mane_summary = raw_dir / "MANE.GRCh38.summary.txt.gz"
    uniprot_mapping = raw_dir / "HUMAN_9606_idmapping.dat.gz"

    missing = []
    if not interpro_xml.exists():
        missing.append(("InterPro XML", interpro_xml, INTERPRO_XML_URL))
    if not mane_summary.exists():
        missing.append(("MANE Summary", mane_summary, MANE_SUMMARY_URL))
    if not uniprot_mapping.exists():
        missing.append(("UniProt Mapping", uniprot_mapping, UNIPROT_MAPPING_URL))

    if missing:
        if download:
            for desc, path, url in missing:
                download_file(url, path, desc)
        else:
            missing_files = [str(p) for _, p, _ in missing]
            raise FileNotFoundError(
                f"Missing required files:\n" + "\n".join(missing_files) +
                "\n\nRun with --download to automatically download them."
            )

    return interpro_xml, mane_summary, uniprot_mapping


@dataclass
class InterproEntry:
    """A single InterPro domain/annotation entry."""
    interpro_id: str
    element_type: str
    match_id: str
    element_name: str
    representative: bool
    positions: list[tuple[int, int]]


@dataclass
class Protein:
    """A protein with its InterPro annotations."""
    uniprot_id: str
    name: str
    length: int
    entries: list[InterproEntry]


def parse_one_match(match_ele: ET.Element) -> Optional[InterproEntry]:
    """Parse a single <match> element into an InterproEntry."""
    grouped = defaultdict(list)
    for ele in match_ele.iter():
        grouped[ele.tag].append(ele)

    assert len(grouped["match"]) == 1
    match_obj = grouped["match"][0]

    # Only keep matches integrated into InterPro (have <ipr> element)
    if "ipr" not in grouped:
        return None

    assert len(grouped["ipr"]) == 1
    ipr = grouped["ipr"][0]

    # Check representative status from <lcn> elements
    is_representative = [x.get("representative") for x in grouped["lcn"]]
    representative = is_representative[0] == "true" if is_representative[0] else False

    # Extract positions
    positions = [(int(x.get("start")), int(x.get("end"))) for x in grouped["lcn"]]

    return InterproEntry(
        interpro_id=ipr.get("id"),
        element_type=ipr.get("type"),
        match_id=match_obj.get("id"),
        element_name=ipr.get("name"),
        representative=representative,
        positions=positions,
    )


def parse_one_protein(prot_ele: ET.Element) -> Protein:
    """Parse a single <protein> element into a Protein object."""
    all_entries = [parse_one_match(x) for x in prot_ele.iter(tag="match")]
    parsed_entries = [x for x in all_entries if x is not None]

    return Protein(
        uniprot_id=prot_ele.get("id"),
        name=prot_ele.get("name"),
        length=int(prot_ele.get("length")),
        entries=parsed_entries,
    )


def parse_interpro_xml(
    input_path: Path,
    target_uniprot_ids: Optional[set[str]] = None,
    print_every: int = 100000,
) -> Generator[Protein, None, None]:
    """
    Stream parse InterPro XML, optionally filtering to target UniProt IDs.

    Args:
        input_path: Path to match_complete.xml.gz
        target_uniprot_ids: If provided, only yield proteins in this set
        print_every: Log progress every N proteins

    Yields:
        Protein objects with their InterPro entries
    """
    logger.info(f"Parsing InterPro XML: {input_path}")
    if target_uniprot_ids:
        logger.info(f"Filtering to {len(target_uniprot_ids)} target UniProt IDs")

    with gzip.open(input_path, "rt") as fh:
        in_block = False
        lines = []
        total_parsed = 0
        total_matched = 0

        for line in fh:
            if line.startswith("<protein "):
                in_block = True
                lines = []

                # Quick check if this protein is in our target set
                if target_uniprot_ids:
                    # Extract ID from line: <protein id="P12345" ...>
                    start = line.find('id="') + 4
                    end = line.find('"', start)
                    protein_id = line[start:end]

                    if protein_id not in target_uniprot_ids:
                        in_block = False
                        continue

            if in_block:
                lines.append(line)

            if line.startswith("</protein>") and in_block:
                in_block = False
                protein = parse_one_protein(ET.fromstringlist(lines))
                total_parsed += 1

                if target_uniprot_ids is None or protein.uniprot_id in target_uniprot_ids:
                    total_matched += 1
                    yield protein

                if total_parsed % print_every == 0:
                    logger.info(f"Processed {total_parsed:,} proteins, matched {total_matched:,}")

        logger.info(f"Finished: {total_parsed:,} proteins processed, {total_matched:,} matched")


def load_mane_summary(mane_path: Path) -> pd.DataFrame:
    """Load MANE summary file."""
    logger.info(f"Loading MANE summary: {mane_path}")
    df = pd.read_csv(mane_path, sep='\t', compression='gzip')
    df = df[df['MANE_status'] == 'MANE Select'].copy()
    logger.info(f"Found {len(df)} MANE Select entries")
    return df


def load_uniprot_mapping(mapping_path: Path, target_ensembl_proteins: set[str]) -> dict[str, str]:
    """
    Load UniProt ID mapping file and create Ensembl protein -> UniProt mapping.

    Args:
        mapping_path: Path to HUMAN_9606_idmapping.dat.gz
        target_ensembl_proteins: Set of Ensembl protein IDs to look up

    Returns:
        Dict mapping Ensembl protein ID -> UniProt accession
    """
    logger.info(f"Loading UniProt ID mapping: {mapping_path}")

    # Build lookup set (include both versioned and unversioned)
    ensp_set = set()
    for ensp in target_ensembl_proteins:
        ensp_set.add(ensp)
        ensp_set.add(ensp.split('.')[0])

    ensp_to_uniprot = {}
    reviewed_accs = set()

    with gzip.open(mapping_path, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue

            uniprot_acc, id_type, id_val = parts

            if id_type == 'Ensembl_PRO':
                id_base = id_val.split('.')[0]
                if id_val in ensp_set or id_base in ensp_set:
                    if id_val not in ensp_to_uniprot:
                        ensp_to_uniprot[id_val] = []
                    ensp_to_uniprot[id_val].append(uniprot_acc)

                    if id_base not in ensp_to_uniprot:
                        ensp_to_uniprot[id_base] = []
                    if uniprot_acc not in ensp_to_uniprot[id_base]:
                        ensp_to_uniprot[id_base].append(uniprot_acc)

            elif id_type == 'UniProtKB-ID' and '_HUMAN' in id_val:
                reviewed_accs.add(uniprot_acc)

    # Select best UniProt accession (prefer reviewed/SwissProt)
    final_mapping = {}
    for ensp, accs in ensp_to_uniprot.items():
        if len(accs) == 1:
            final_mapping[ensp] = accs[0]
        else:
            reviewed = [a for a in accs if a in reviewed_accs]
            final_mapping[ensp] = reviewed[0] if reviewed else accs[0]

    logger.info(f"Mapped {len(final_mapping)} Ensembl proteins to UniProt")
    return final_mapping


def extract_domains(
    protein: Protein,
    use_representative: bool = True,
    exclude_types: set[str] = EXCLUDED_TYPES,
) -> list[dict]:
    """
    Extract domain records from a protein.

    Args:
        protein: Protein object with InterPro entries
        use_representative: Only use representative domains for Domain/Family/Repeat
        exclude_types: Domain types to exclude

    Returns:
        List of domain record dicts
    """
    records = []

    for entry in protein.entries:
        # Skip excluded types
        if entry.element_type in exclude_types:
            continue

        # For Domain/Family/Repeat, only use representative if requested
        if use_representative and entry.element_type in REPRESENTATIVE_TYPES:
            if not entry.representative:
                continue

        for start, end in entry.positions:
            records.append({
                'protein_id_uniprot': protein.uniprot_id,
                'domain_id_interpro': entry.interpro_id,
                'domain_name': entry.element_name,
                'domain_type': entry.element_type,
                'source_db': entry.match_id,
                'domain_start_aa': start,
                'domain_end_aa': end,
                'representative': entry.representative,
            })

    return records


def main():
    parser = argparse.ArgumentParser(
        description='Parse InterPro domains for MANE Select transcripts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # With automatic downloads (recommended for first run):
  python parse_interpro_domains.py --data-dir ./data --output mane_domains.parquet --download

  # With existing files:
  python parse_interpro_domains.py \\
      --interpro-xml /path/to/match_complete.xml.gz \\
      --mane-summary /path/to/MANE.GRCh38.summary.txt.gz \\
      --uniprot-mapping /path/to/HUMAN_9606_idmapping.dat.gz \\
      --output mane_domains.parquet
        """
    )

    # Data source options (either --data-dir or individual paths)
    data_group = parser.add_argument_group('Data Sources')
    data_group.add_argument(
        '--data-dir',
        type=Path,
        help='Directory for data files (will look for/download files in data-dir/raw/)'
    )
    data_group.add_argument(
        '--download',
        action='store_true',
        help='Download missing data files automatically'
    )
    data_group.add_argument(
        '--interpro-xml',
        type=Path,
        help='Path to match_complete.xml.gz from InterPro'
    )
    data_group.add_argument(
        '--mane-summary',
        type=Path,
        help='Path to MANE.GRCh38.summary.txt.gz'
    )
    data_group.add_argument(
        '--uniprot-mapping',
        type=Path,
        help='Path to HUMAN_9606_idmapping.dat.gz'
    )

    # Output options
    output_group = parser.add_argument_group('Output')
    output_group.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output parquet file path'
    )

    # Processing options
    proc_group = parser.add_argument_group('Processing Options')
    proc_group.add_argument(
        '--include-non-representative',
        action='store_true',
        help='Include non-representative domains (default: only representative for Domain/Family/Repeat)'
    )

    args = parser.parse_args()

    # Resolve data file paths
    if args.data_dir:
        interpro_xml, mane_summary, uniprot_mapping = ensure_data_files(
            args.data_dir, download=args.download
        )
    elif args.interpro_xml and args.mane_summary and args.uniprot_mapping:
        interpro_xml = args.interpro_xml
        mane_summary = args.mane_summary
        uniprot_mapping = args.uniprot_mapping
    else:
        parser.error("Either --data-dir or all three file paths (--interpro-xml, --mane-summary, --uniprot-mapping) are required")

    # Load MANE summary
    mane_df = load_mane_summary(mane_summary)

    # Get Ensembl protein IDs from MANE
    ensembl_proteins = set(mane_df['Ensembl_prot'].dropna().tolist())
    logger.info(f"Found {len(ensembl_proteins)} Ensembl protein IDs in MANE")

    # Load UniProt mapping
    ensp_to_uniprot = load_uniprot_mapping(uniprot_mapping, ensembl_proteins)

    # Create reverse mapping: UniProt -> Ensembl transcript
    uniprot_to_transcript = {}
    for _, row in mane_df.iterrows():
        ensp = row.get('Ensembl_prot')
        enst = row.get('Ensembl_nuc')
        gene = row.get('symbol')

        if pd.notna(ensp):
            uniprot = ensp_to_uniprot.get(ensp) or ensp_to_uniprot.get(ensp.split('.')[0])
            if uniprot:
                uniprot_to_transcript[uniprot] = {
                    'transcript_id_ensembl': enst,
                    'protein_id_ensembl': ensp,
                    'gene_symbol': gene,
                }

    logger.info(f"Created mapping for {len(uniprot_to_transcript)} UniProt -> transcript pairs")

    # Target UniProt IDs to extract
    target_uniprot_ids = set(uniprot_to_transcript.keys())

    # Parse InterPro XML and extract domains
    all_records = []
    proteins_with_domains = 0

    for protein in parse_interpro_xml(interpro_xml, target_uniprot_ids):
        domains = extract_domains(
            protein,
            use_representative=not args.include_non_representative
        )

        if domains:
            proteins_with_domains += 1
            # Add transcript info to each domain record
            transcript_info = uniprot_to_transcript.get(protein.uniprot_id, {})
            for record in domains:
                record.update(transcript_info)
            all_records.extend(domains)

    logger.info(f"Extracted {len(all_records)} domain records from {proteins_with_domains} proteins")

    # Create DataFrame and save
    if all_records:
        df = pd.DataFrame(all_records)

        # Reorder columns
        column_order = [
            'transcript_id_ensembl',
            'protein_id_ensembl',
            'protein_id_uniprot',
            'gene_symbol',
            'domain_id_interpro',
            'domain_name',
            'domain_type',
            'source_db',
            'domain_start_aa',
            'domain_end_aa',
            'representative',
        ]
        df = df[[c for c in column_order if c in df.columns]]

        # Sort
        df = df.sort_values(['transcript_id_ensembl', 'domain_start_aa'])

        # Save
        df.to_parquet(args.output, index=False, compression='snappy')
        logger.info(f"Saved {len(df)} records to {args.output}")

        # Print summary
        logger.info("\n=== Summary ===")
        logger.info(f"Total domain records: {len(df)}")
        logger.info(f"Unique proteins: {df['protein_id_uniprot'].nunique()}")
        logger.info(f"Unique InterPro entries: {df['domain_id_interpro'].nunique()}")
        logger.info(f"Unique genes: {df['gene_symbol'].nunique()}")
        logger.info(f"\nDomain types:")
        for dtype, count in df['domain_type'].value_counts().items():
            logger.info(f"  {dtype}: {count}")
    else:
        logger.warning("No domain records extracted!")


if __name__ == '__main__':
    main()
