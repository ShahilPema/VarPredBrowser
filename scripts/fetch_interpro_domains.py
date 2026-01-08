#!/usr/bin/env python3
"""
Fetch InterPro representative domains from all member databases.

The /entry/all/ endpoint doesn't include all databases (missing ssf, smart, panther).
This script fetches from both /entry/all/ and the individual database endpoints
to get complete representative domain coverage.
"""

import json
import requests
import time
from pathlib import Path
from tqdm import tqdm

# Databases not included in /entry/all/ that need separate queries
EXTRA_DATABASES = ['ssf', 'smart', 'panther']

# Only keep these entry types (exclude 'family', 'repeat', etc.)
VALID_DOMAIN_TYPES = {'domain', 'homologous_superfamily'}


def fetch_from_endpoint(url, delay=0.1):
    """Fetch domains from a single API endpoint."""
    try:
        response = requests.get(url, headers={'Accept': 'application/json'}, timeout=30)

        if response.status_code == 200:
            time.sleep(delay)
            return response.json()
        elif response.status_code == 404:
            time.sleep(delay)
            return {'results': []}
        else:
            time.sleep(delay)
            return None

    except requests.exceptions.RequestException:
        return None


def parse_domains(data, only_representative=True):
    """Parse domains from API response."""
    domains = []

    for result in data.get('results', []):
        metadata = result.get('metadata', {})
        acc_id = metadata.get('accession', '')
        name = metadata.get('name', '') or ''
        entry_type = metadata.get('type', '')
        source_db = metadata.get('source_database', '')

        # Skip entries that aren't domains or homologous superfamilies
        if entry_type not in VALID_DOMAIN_TYPES:
            continue

        # Get InterPro ID if this signature is integrated
        integrated = metadata.get('integrated', None)
        ipr_id = integrated if integrated else ''

        for protein in result.get('proteins', []):
            for loc in protein.get('entry_protein_locations', []):
                is_rep = loc.get('representative', False)

                if not only_representative or is_rep:
                    for frag in loc.get('fragments', []):
                        domains.append({
                            'member_db_id': acc_id,
                            'interpro_id': ipr_id,
                            'domain_name': name,
                            'domain_type': entry_type,
                            'source_db': source_db,
                            'start_aa': frag.get('start', 0),
                            'end_aa': frag.get('end', 0),
                            'representative': is_rep
                        })

    return domains


def fetch_representative_domains(acc, delay=0.1):
    """Fetch representative domains for a single protein from all member databases."""
    all_domains = []

    # First fetch from /entry/all/ (includes cathgene3d, cdd, pfam, prints, profile, interpro)
    url_all = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}"
    data = fetch_from_endpoint(url_all, delay=delay)

    if data is None:
        return None  # Request failed

    all_domains.extend(parse_domains(data, only_representative=True))

    # Then fetch from extra databases not included in /entry/all/
    for db in EXTRA_DATABASES:
        url_db = f"https://www.ebi.ac.uk/interpro/api/entry/{db}/protein/uniprot/{acc}"
        data = fetch_from_endpoint(url_db, delay=delay)

        if data is not None:
            all_domains.extend(parse_domains(data, only_representative=True))

    return all_domains


def main():
    # Load existing cache to get protein list
    cache_path = Path('data/cache/interpro_domains_cache.json')
    rep_cache_path = Path('data/cache/interpro_representative_domains.json')

    # Create cache directory if needed
    rep_cache_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading existing cache for protein list...")
    with open(cache_path, 'r') as f:
        old_cache = json.load(f)

    proteins = list(old_cache.keys())
    print(f"Will fetch representative domains for {len(proteins)} proteins...")
    print(f"Fetching from /entry/all/ + extra databases: {EXTRA_DATABASES}")

    new_cache = {}
    failed = []

    for acc in tqdm(proteins, desc="Fetching representative domains"):
        domains = fetch_representative_domains(acc)

        if domains is None:
            failed.append(acc)
        elif domains:  # Only store if there are representative domains
            new_cache[acc] = domains

    # Save new cache
    print(f"\nSaving representative domains cache to {rep_cache_path}")
    with open(rep_cache_path, 'w') as f:
        json.dump(new_cache, f)

    print(f"\nDone! Found representative domains for {len(new_cache)} proteins.")
    if failed:
        print(f"Failed to fetch {len(failed)} proteins")

if __name__ == '__main__':
    main()
