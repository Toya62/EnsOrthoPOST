#!/usr/bin/env python3
"""
Extract ortholog clusters from Ensembl for target species
"""

import requests
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote


def fetch_single_gene_orthologs(human_id, symbol, throttle_delay, allowed_species=None):
    """
    Fetch orthologs for a single human gene via Ensembl API.
    
    Args:
        human_id: Ensembl human gene ID
        symbol: Human gene symbol
        throttle_delay: Seconds to wait before request (rate limiting)
    
    Returns:
        tuple: (human_id, symbol, orthologs_dict) or (human_id, symbol, None, error_msg) on failure
    """
    encoded_symbol = quote(symbol, safe='')
    urls = [
        f"https://rest.ensembl.org/homology/symbol/homo_sapiens/{encoded_symbol}?content-type=application/json",
        f"https://rest.ensembl.org/homology/id/homo_sapiens/{human_id}?content-type=application/json",
    ]
    headers = {"Content-Type": "application/json"}

    max_retries = 3
    last_error = None
    for url in urls:
        for attempt in range(max_retries):
            # Apply throttle delay before request
            time.sleep(throttle_delay)

            try:
                response = requests.get(url, headers=headers, timeout=30)

                # Handle rate limiting
                if response.status_code == 429:
                    wait_time = float(response.headers.get("Retry-After", 2))
                    print(f"  Rate limited, waiting {wait_time}s...", flush=True)
                    time.sleep(wait_time)
                    continue

                # Retry on server errors (5xx)
                if 500 <= response.status_code < 600:
                    if attempt < max_retries - 1:
                        wait_time = 2 ** attempt  # exponential backoff
                        print(f"  Server error {response.status_code}, retrying in {wait_time}s...", flush=True)
                        time.sleep(wait_time)
                        continue
                    last_error = f"{response.status_code} Server Error (max retries exceeded)"
                    break

                # Symbol lookups occasionally fail on aliases; fallback to id endpoint
                if response.status_code in (400, 404):
                    last_error = f"{response.status_code} {response.reason}"
                    break

                # Other errors - fail immediately
                response.raise_for_status()
                data = response.json()

                # Initialize the cluster for this gene
                gene_orthos = defaultdict(list)

                # Extract homologies from response
                if 'data' in data and data['data']:
                    homologies = data['data'][0].get('homologies', [])

                    # Filter to orthologues (not paralogues)
                    for entry in homologies:
                        entry_type = entry.get('type', '')

                        # Only keep one2one and one2many ortholog types
                        if entry_type not in ['ortholog_one2one', 'ortholog_one2many']:
                            continue

                        species = entry.get('target', {}).get('species')
                        ortho_id = entry.get('target', {}).get('id')

                        if allowed_species is not None and species not in allowed_species:
                            continue

                        if species and ortho_id:
                            gene_orthos[species].append(ortho_id)

                return (human_id, symbol, dict(gene_orthos))

            except Exception as e:
                last_error = str(e)
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    time.sleep(wait_time)
                    continue
                break
    
    return (human_id, symbol, None, last_error or "Max retries exceeded")


def extract_ortho_clusters(human_gene_map, target_mammals, max_workers, throttle_delay):
    """
    Fetch orthologs for human genes via Ensembl API using parallel workers.
    
    Args:
        human_gene_map: dict {ensembl_id: symbol}
        target_mammals: list of target species (for validation)
        max_workers: number of parallel workers
        throttle_delay: seconds to wait between API requests
    
    Returns:
        dict: {human_id: {'symbol': str, 'orthologs': {species: [ids]}}}
    """
    ortho_clusters = {}
    total_genes = len(human_gene_map)
    processed = 0
    failed_genes = []
    ortholog_issues = []
    
    print(f"Extracting orthologs for {total_genes} genes using {max_workers} workers (throttle: {throttle_delay}s)...", flush=True)

    allowed_species = set(target_mammals) if target_mammals else None

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs to the executor
        futures = {
            executor.submit(fetch_single_gene_orthologs, human_id, symbol, throttle_delay, allowed_species): 
            (human_id, symbol) 
            for human_id, symbol in human_gene_map.items()
        }
        
        # Process results as they complete
        for future in as_completed(futures):
            processed += 1
            
            # Print progress every 1000 genes
            if processed % 1000 == 0 or processed == 1:
                print(f"Progress: {processed}/{total_genes} genes processed...", flush=True)
            
            try:
                result = future.result()
                human_id, symbol, orthologs = result[0], result[1], result[2]
                
                if orthologs is not None:
                    # Add the human gene itself to the ortholog list
                    if 'homo_sapiens' not in orthologs:
                        orthologs['homo_sapiens'] = []
                    orthologs['homo_sapiens'].append(human_id)

                    non_human_ortholog_count = sum(
                        len(ids) for species, ids in orthologs.items() if species != 'homo_sapiens'
                    )
                    if non_human_ortholog_count == 0:
                        ortholog_issues.append({
                            'ensembl_id': human_id,
                            'symbol': symbol,
                            'issue_type': 'no_nonhuman_orthologs',
                            'reason': 'No ortholog_one2one/ortholog_one2many hits found in target species list'
                        })
                    
                    ortho_clusters[human_id] = {
                        'symbol': symbol,
                        'orthologs': orthologs
                    }
                else:
                    error_msg = result[3] if len(result) > 3 else "Unknown error"
                    print(f"   Error fetching {symbol} ({human_id}): {error_msg}", flush=True)
                    failed_genes.append({
                        'ensembl_id': human_id,
                        'symbol': symbol,
                        'reason': error_msg,
                    })
                    
            except Exception as e:
                human_id, symbol = futures[future]
                print(f"   Exception for {symbol} ({human_id}): {e}", flush=True)
                failed_genes.append({
                    'ensembl_id': human_id,
                    'symbol': symbol,
                    'reason': str(e),
                })
    
    print(f"\n Completed cluster extraction for {len(ortho_clusters)} human genes.", flush=True)
    if failed_genes:
        print(f" Failed to fetch orthologs for {len(failed_genes)} genes", flush=True)
    if ortholog_issues:
        print(f" {len(ortholog_issues)} genes have no non-human orthologs in target species", flush=True)
    
    # Store failed genes for statistics
    ortho_clusters['_failed_genes'] = failed_genes
    ortho_clusters['_ortholog_issues'] = ortholog_issues
    
    return ortho_clusters
