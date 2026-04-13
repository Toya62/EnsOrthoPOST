#!/usr/bin/env python3
"""
Extract human genes from HGNC for a gene group
"""

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import csv
from io import StringIO


def extract_control_genes(group_id):
    """
    Fetch protein-coding genes from HGNC for specified gene group ID.
    
    Args:
        group_id: HGNC gene group ID (e.g., 141, 149, 159, 595)
    
    Returns:
        dict: {ensembl_id: symbol}
    """
    
    print(f"Extracting genes from HGNC gene group ID {group_id}...")
    
    # Setup: headers and session with retries
    session = requests.Session()
    retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503], allowed_methods=("GET",))
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    
    # Try the download endpoint (HGNC gene group download interface)
    genes = {}
    
    print(f"  Attempting download endpoint...")
    try:
        url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch&format=json"
        response = session.get(url, timeout=15)
        response.raise_for_status()
        
        if response.text:
            content_type = (response.headers.get("content-type") or "").lower()

            # Preferred path: JSON format from HGNC download endpoint
            if "application/json" in content_type or response.text.lstrip().startswith("["):
                rows = response.json()
                if isinstance(rows, list):
                    for row in rows:
                        ensembl_id = (row.get("ensemblGeneID") or "").strip()
                        symbol = (row.get("approvedSymbol") or "").strip()
                        locus_type = (row.get("locusType") or "").lower()

                        # Only keep protein-coding genes
                        if ensembl_id and symbol and "protein" in locus_type:
                            genes[ensembl_id] = symbol
            else:
                # Legacy fallback: parse tab-delimited text download format
                reader = csv.DictReader(StringIO(response.text), delimiter='\t')
                for row in reader:
                    ensembl_id = row.get('Ensembl gene ID', '').strip() if row.get('Ensembl gene ID') else ''
                    symbol = row.get('Approved symbol', '').strip() if row.get('Approved symbol') else ''
                    locus_type = (row.get('Locus type', '') or '').lower()

                    # Only keep protein-coding genes
                    if ensembl_id and symbol and 'protein' in locus_type:
                        genes[ensembl_id] = symbol

            if genes:
                print(f"OK: Found {len(genes)} protein-coding genes")
                return genes
    except Exception as e:
        print(f"WARNING: Download endpoint failed: {e}")
    
    # Fallback to REST API Search
    print(f"  Attempting REST API search...")
    headers = {"Accept": "application/json"}
    try:
        hgnc_url = f"https://rest.genenames.org/search/gene_group_id:{group_id}"
        hgnc_resp = session.get(hgnc_url, headers=headers, timeout=15)
        if hgnc_resp.status_code == 200:
            hgnc_data = hgnc_resp.json()
            docs = hgnc_data.get('response', {}).get('docs', [])
            
            if docs:
                symbols = [doc.get('symbol') for doc in docs if doc.get('symbol')]
                
                # Look up each symbol in Ensembl to get IDs and verify protein-coding
                for symbol in symbols:
                    try:
                        ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?expand=0"
                        ens_resp = session.get(ensembl_url, headers=headers, timeout=10)
                        if ens_resp.status_code == 200:
                            ens_data = ens_resp.json()
                            ensembl_id = ens_data.get('id')
                            locus_type = (ens_data.get('biotype', '') or '').lower()
                            if ensembl_id and 'protein' in locus_type:
                                genes[ensembl_id] = symbol
                    except:
                        pass
                
                if genes:
                    print(f"OK: Found {len(genes)} protein-coding genes")
                    return genes
    except Exception as e:
        print(f"WARNING: REST API search failed: {e}")
    
        # No genes found
        print(f"WARNING: No genes found for group ID {group_id}")
    print(f"           (Group {group_id} may not be indexed in HGNC or may have no protein-coding genes)")
    return {}
