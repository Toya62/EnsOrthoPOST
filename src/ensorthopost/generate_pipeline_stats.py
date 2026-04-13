#!/usr/bin/env python3
"""
Pipeline Statistics Generator for GLUH EnsOrthoPOST
Analyzes pipeline output and generates comprehensive statistics for research.
"""
import json
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd


def load_fasta_cluster_map(output_dir, group_id):
    """Load cluster assignments from FASTA headers.

    Returns:
        dict: {(species, member_ensembl_id): cluster_id}
    """
    base = Path(output_dir)
    fasta_file = base / "fasta" / f"{group_id}.fa"
    if not fasta_file.exists():
        return {}

    cluster_map = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            if 'cluster=' not in header or 'member=' not in header or 'species=' not in header:
                continue

            parts = header.split(';')
            values = {}
            for part in parts:
                if '=' in part:
                    key, value = part.split('=', 1)
                    values[key] = value

            cluster_id = values.get('cluster')
            member_id = values.get('member')
            species = values.get('species')
            if cluster_id and member_id and species:
                cluster_map[(species, member_id)] = cluster_id

    return cluster_map


def load_metadata_file(output_dir, group_id):
    """Load ortholog metadata JSON file"""
    base = Path(output_dir)
    metadata_file = base / "stats" / f"{group_id}_ortholog_metadata.json"
    if not metadata_file.exists():
        metadata_file = base / "stats" / f"ortholog_metadata_{group_id}.json"
    if not metadata_file.exists():
        metadata_file = base / f"{group_id}_ortholog_metadata.json"
    if not metadata_file.exists():
        metadata_file = base / f"ortholog_metadata_{group_id}.json"
    if not metadata_file.exists():
        print(f"WARNING: Metadata file not found: {metadata_file}")
        return []
    
    with open(metadata_file, 'r') as f:
        return json.load(f)


def _extract_transcript_id(attr_str):
    """Extract transcript_id from a GTF attributes field."""
    for attr in attr_str.split(';'):
        attr = attr.strip()
        if 'transcript_id' not in attr:
            continue

        if '"' in attr:
            parts = attr.split('"')
            if len(parts) >= 2:
                return parts[1]
        elif '=' in attr:
            return attr.split('=', 1)[1].strip()

    return None


def create_gene_stats_dataframe(output_dir, group_id):
    """
    Create pandas DataFrame with per-gene statistics.
    
    Columns:
        - group_id: Gene group ID
        - cluster: Source human Ensembl gene ID from FASTA header
        - human_gene_id: Human Ensembl gene ID (cluster source)
        - human_gene_symbol: Human gene symbol
        - ensembl_id_gene: Ortholog Ensembl gene ID
        - ensembl_id_transcript: Ensembl transcript ID (from GTF)
        - gene_length_bp: Gene length in base pairs
        - species: Species name
    
    Returns:
        pd.DataFrame or None if no data available
    """
    
    metadata = load_metadata_file(output_dir, group_id)
    if not metadata:
        print("WARNING: No metadata available for dataframe generation")
        return None
    
    # Load ortholog clusters to map orthologs back to human genes
    base = Path(output_dir)
    fasta_cluster_map = load_fasta_cluster_map(output_dir, group_id)
    cluster_file = base / "stats" / f"{group_id}_ortho_clusters.json"
    if not cluster_file.exists():
        cluster_file = base / "stats" / f"ortho_clusters_{group_id}.json"
    if not cluster_file.exists():
        cluster_file = base / f"{group_id}_ortho_clusters.json"
    if not cluster_file.exists():
        cluster_file = base / f"ortho_clusters_{group_id}.json"
    ortholog_to_human = {}  # {(species, ortholog_id): (human_gene_id, human_symbol)}
    
    if cluster_file.exists():
        with open(cluster_file, 'r') as f:
            ortho_clusters = json.load(f)
        
        for human_gene_id, gene_data in ortho_clusters.items():
            if human_gene_id.startswith('_') or not isinstance(gene_data, dict):
                continue
            
            human_symbol = gene_data.get('symbol', '')
            orthologs = gene_data.get('orthologs', {})
            
            for species, ortho_ids in orthologs.items():
                for ortho_id in ortho_ids:
                    key = (species, ortho_id)
                    if key not in ortholog_to_human:
                        ortholog_to_human[key] = (human_gene_id, human_symbol)
    else:
        print("WARNING: Ortholog cluster file not found - human gene mapping unavailable")

    # Fallback: if FASTA mapping is unavailable for a row, use the human gene ID.
    cluster_id_map = {}
    if cluster_file.exists():
        human_gene_ids = [k for k in ortho_clusters.keys() if not str(k).startswith('_')]
        for human_gene_id in human_gene_ids:
            cluster_id_map[human_gene_id] = str(human_gene_id)
    
    # Parse GTF to extract transcript IDs
    gtf_file = base / "gtf" / f"{group_id}.gtf"
    if not gtf_file.exists():
        gtf_file = base / f"{group_id}.gtf"
    transcript_map = defaultdict(set)  # {ensembl_id: {transcript_ids}}
    
    if gtf_file.exists():
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Extract transcript_id from attributes (9th field)
                attr_str = fields[8]
                transcript_id = _extract_transcript_id(attr_str)
                
                # Extract gene ID from first field
                seq_id = fields[0]
                if 'member=' in seq_id:
                    ensembl_id = seq_id.split('member=', 1)[1].split(';', 1)[0]
                elif '_' in seq_id:
                    ensembl_id = '_'.join(seq_id.split('_')[1:])
                else:
                    ensembl_id = seq_id
                
                # Always map - use gene_id as fallback if no transcript_id
                if not transcript_id:
                    transcript_id = ensembl_id
                
                transcript_map[ensembl_id].add(transcript_id)
    
    # Build dataframe rows
    rows = []
    for meta in metadata:
        ensembl_id = meta['ensembl_id']
        species = meta['species']
        length = meta.get('length', 0)
        
        # Get human gene info for this ortholog
        human_info = ortholog_to_human.get((species, ensembl_id), (None, None))
        human_gene_id, human_symbol = human_info

        # Get cluster from FASTA header (authoritative ETL output); fallback to derived map
        cluster_id = fasta_cluster_map.get((species, ensembl_id), cluster_id_map.get(human_gene_id, ""))
        
        # Get transcript IDs for this gene (always non-empty now)
        transcripts = transcript_map.get(ensembl_id, {ensembl_id})
        
        for transcript_id in transcripts:
            rows.append({
                'group_id': group_id,
                'cluster': cluster_id,
                'human_gene_id': human_gene_id if human_gene_id else '',
                'human_gene_symbol': human_symbol if human_symbol else '',
                'ensembl_id_gene': ensembl_id,
                'ensembl_id_transcript': transcript_id if transcript_id else ensembl_id,  # Fallback to gene ID
                'gene_length_bp': int(length),
                'species': species
            })
    
    if not rows:
        print("WARNING: No genes found for dataframe")
        return None
    
    df = pd.DataFrame(rows)
    return df


def save_gene_stats_dataframe(df, output_dir, group_id):
    """Save gene statistics dataframe to CSV"""
    if df is None:
        print("WARNING: No dataframe to save")
        return
    # Save CSV
    stats_dir = Path(output_dir) / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    csv_file = stats_dir / f"{group_id}_gene_stats.csv"
    df.to_csv(csv_file, index=False)
    print(f"OK: Gene statistics saved: {csv_file}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate simple consolidated statistics for GLUH EnsOrthoPOST outputs'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Pipeline output directory'
    )
    parser.add_argument(
        '--group-id',
        type=str,
        required=True,
        help='Gene group ID (OR, TAAR, etc.)'
    )
    args = parser.parse_args()
    # Generate and save gene statistics dataframe (CSV output only)
    print("\nGenerating gene statistics dataframe...")
    gene_df = create_gene_stats_dataframe(args.output_dir, args.group_id)
    if gene_df is not None:
        save_gene_stats_dataframe(gene_df, args.output_dir, args.group_id)
    else:
        print("WARNING: Could not generate gene statistics dataframe")
    return 0

if __name__ == '__main__':
    sys.exit(main())
