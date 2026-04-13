#!/usr/bin/env python3
"""
EnsOrthoPOST Pipeline - Extract, Transform, Load ortholog sequences
"""

import sys
import json
from pathlib import Path

from .extract_control_genes import extract_control_genes
from .helper_define_species import load_species_list
from .extract_ortho_clusters import extract_ortho_clusters
from .gene_metadata_report import generate_gene_metadata
from .build_corpus_files import build_corpus_files
from .generate_pipeline_stats import (
    create_gene_stats_dataframe, save_gene_stats_dataframe
)


def run_main_orchestrator(
    group_id,
    output_dir=None,
    max_workers=2,
    throttle=0.3,
    species_file=None,
    include_human=False,
    deterministic=True,
):
    """
    EnsOrthoPOST Pipeline - Extract, Transform, Load ortholog sequences
    
    Args:
        group_id: HGNC numeric gene-group ID (e.g., 141, 159)
        output_dir: output directory (optional)
        max_workers: parallel workers (default: 2)
        throttle: API throttle in seconds (default: 0.3)
        species_file: path to species file (.txt)
        include_human: include human entries in per-family FASTA/GTF outputs
    """
    
    # Setup output directory structure
    if output_dir is None:
        output_dir = Path("/work/welsch/raw_data/families_test")
    else:
        output_dir = Path(output_dir)

    output_dir = output_dir
    (output_dir / "fasta").mkdir(parents=True, exist_ok=True)
    (output_dir / "gtf").mkdir(parents=True, exist_ok=True)
    (output_dir / "stats").mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}\nEnsOrthoPOST Pipeline: Group {group_id} | Workers={max_workers}\n{'='*60}\n")
    
    # STEP 1: Extract human genes
    print("[1/6] Extracting human genes from HGNC...")
    human_genes = extract_control_genes(group_id)
    if not human_genes:
        print(f"ERROR: No human genes found for group {group_id} - invalid ID or HGNC API unavailable.")
        # Write failure report
        (output_dir / "stats").mkdir(parents=True, exist_ok=True)
        failure_report = {
            'group_id': group_id,
            'reason': 'No human genes found (invalid ID or HGNC API unavailable)'
        }
        failure_file = output_dir / "stats" / f"{group_id}_failed_genes.json"
        with open(failure_file, 'w') as f:
            json.dump(failure_report, f, indent=2)
        print(f"INFO: Failure report saved to: {failure_file}\n")
        return 1
    print(f"OK: {len(human_genes)} genes\n")

    # Save input human anchor genes (HGNC group members) only; excludes orthologs
    genes_list_file = output_dir / "genes_list.txt"
    with open(genes_list_file, "w") as f:
        f.write("\n".join(
            f"{ensembl_id}\t{symbol}" for ensembl_id, symbol in human_genes.items()
        ))
    print(f"INFO: Genes list saved to: {genes_list_file.name}\n")

    # STEP 2: Load species
    print("[2/6] Loading target species...")
    try:
        species_list = load_species_list(species_file)
    except Exception as e:
        print(f"ERROR: Species selection failed: {e}")
        return 1

    print(f"  - Input species: {len(set(species_list))}")
    print(f"  - Selected species: {len(species_list)}")

    print(f"OK: {len(species_list)} species\n")

    # Save species list to output folder
    species_list_file = output_dir / "species_list.txt"
    with open(species_list_file, "w") as f:
        f.write("\n".join(species_list))
    print(f"INFO: Species list saved to: {species_list_file.name}\n")

    # STEP 3: Fetch orthologs
    print("[3/6] Fetching orthologs...")
    ortho_clusters = extract_ortho_clusters(human_genes, species_list, max_workers, throttle)

    # Check for real ortholog clusters (exclude internal keys)
    num_real = len([k for k in ortho_clusters.keys() if not str(k).startswith('_')])
    if num_real == 0:
        print(f"WARNING: No ortholog clusters found for group {group_id}: skipping (single-species or no orthologs).")
        # Write failure report
        (output_dir / "stats").mkdir(parents=True, exist_ok=True)
        failure_report = {
            'group_id': group_id,
            'reason': 'No ortholog clusters found (single-species or no orthologs)'
        }
        failure_file = output_dir / "stats" / f"{group_id}_failed_genes.json"
        with open(failure_file, 'w') as f:
            json.dump(failure_report, f, indent=2)
        print(f"INFO: Failure report saved to: {failure_file}\n")
        return 2   # return 2 so sbatch can log as SKIPPED

    # Internal keys start with '_' and are not real biological clusters
    num_clusters = num_real
    print(f"{num_clusters} ortholog clusters (one per human gene; humans are included in each cluster)\n")
    
    # (Removed unreachable dead code: if not ortho_clusters)
    
    # STEP 4: Get metadata
    print("[4/6] Fetching metadata...")
    metadata = generate_gene_metadata(ortho_clusters, output_dir, group_id, max_workers, throttle)
    print(f"OK: Metadata saved\n")
    
    # STEP 5: Build corpus
    print("[5/6] Building FASTA + GTF...")
    build_corpus_files(
        ortho_clusters,
        metadata,
        output_dir,
        group_id,
        max_workers,
        deterministic=deterministic,
        include_human=include_human,
    )
    print(f"OK: Done\n")
    
    # Save ortholog clusters for downstream analysis
    cluster_file = output_dir / "stats" / f"{group_id}_ortho_clusters.json"
    ortho_clusters_clean = {
        k: v for k, v in ortho_clusters.items() if not str(k).startswith('_')
    }
    with open(cluster_file, 'w') as f:
        json.dump(ortho_clusters_clean, f, indent=2)
    print(f"OK: Saved ortholog clusters to {cluster_file.name}\n")
    
    # Save failure report for downstream stats
    failed_genes = ortho_clusters.get('_failed_genes', [])
    ortholog_issues = ortho_clusters.get('_ortholog_issues', [])
    if failed_genes or ortholog_issues:
        normalized_failed = []
        for item in failed_genes:
            if isinstance(item, dict):
                normalized_failed.append({
                    'ensembl_id': item.get('ensembl_id', ''),
                    'symbol': item.get('symbol', ''),
                    'reason': item.get('reason', 'Unknown error'),
                })
            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                normalized_failed.append({
                    'ensembl_id': item[0],
                    'symbol': item[1],
                    'reason': 'Unknown error',
                })

        failure_report = {
            'failed_count': len(normalized_failed),
            'failed_genes': normalized_failed,
            'missing_nonhuman_ortholog_count': len(ortholog_issues),
            'missing_nonhuman_ortholog_genes': ortholog_issues,
        }
        failure_file = output_dir / "stats" / f"{group_id}_failed_genes.json"
        with open(failure_file, 'w') as f:
            json.dump(failure_report, f, indent=2)
        print(f"INFO: Failed/issue report saved to: {failure_file}\n")
    
    # STEP 6: Generate gene statistics CSV (even if some API failures occurred)
    print("[6/6] Generating gene statistics CSV...")
    try:
        gene_df = create_gene_stats_dataframe(str(output_dir), group_id)
        if gene_df is not None:
            save_gene_stats_dataframe(gene_df, str(output_dir), group_id)
        print("OK: Gene statistics complete\n")
    except Exception as e:
        print(f"ERROR: Gene stats generation failed: {e}\n")

    print(f"\n{'='*60}\nDONE: {output_dir}\n{'='*60}\n")
    return 0


def main():
    import argparse

    parser = argparse.ArgumentParser(description='EnsOrthoPOST Pipeline - Ortholog Extraction')
    parser.add_argument('--group-id', type=int, default=141, help='Numeric group ID (e.g., 141, 215, etc)')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory (required)')
    parser.add_argument('--workers', type=int, default=2, help='Parallel workers (2-4)')
    parser.add_argument('--throttle', type=float, default=0.3, help='API throttle (sec)')
    parser.add_argument('--species-file', type=str, required=True, help='Path to species list (.txt), one species per line')
    parser.add_argument('--include-human', dest='include_human', action='store_true', default=False,
                        help='Include human anchor genes in per-family FASTA and GTF output (default: excluded)')
    parser.add_argument('--deterministic', dest='deterministic', action='store_true', default=True,
                        help='Use deterministic per-gene flanking buffer sizes (default)')
    parser.add_argument('--non-deterministic', dest='deterministic', action='store_false',
                        help='Use randomized flanking buffer sizes')

    args = parser.parse_args()

    return run_main_orchestrator(
        args.group_id,
        args.output_dir,
        args.workers,
        args.throttle,
        args.species_file,
        args.include_human,
        args.deterministic,
    )


if __name__ == '__main__':
    sys.exit(main())
