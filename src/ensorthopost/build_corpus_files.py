#!/usr/bin/env python3
"""
Build FASTA and GTF corpus files with buffered sequences
"""

import re
import sys
import hashlib
import threading
import requests
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Import Tiberius reformat functions
sys.path.insert(0, str(Path(__file__).parent.parent / 'Tiberius' / 'tiberius'))
from reformat_gtf import parse_attributes_gff3, format_attributes_gtf


_THREAD_LOCAL = threading.local()


def _create_session():
    session = requests.Session()
    retry = Retry(
        total=3,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503],
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def _get_thread_session():
    if not hasattr(_THREAD_LOCAL, "session"):
        _THREAD_LOCAL.session = _create_session()
    return _THREAD_LOCAL.session


def _format_fasta_sequence(seq, width=80):
    """Wrap FASTA sequence lines to a fixed width."""
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def build_corpus_files(
    ortho_clusters,
    metadata,
    output_dir,
    group_id,
    max_workers,
    deterministic=True,
    include_human=False,
):
    """
    Build FASTA sequences and GTF feature rows with flanking buffers.

    Args:
        ortho_clusters: dictionary of ortholog clusters
        metadata: list of metadata records
        output_dir: Path to output directory
        group_id: gene group ID
        max_workers: number of parallel workers
        deterministic: use deterministic per-gene flanking values
        include_human: include human entries in per-family FASTA/GTF outputs
    """

    # Create metadata lookup
    metadata_lookup = {(m['ensembl_id'], m['species']): m for m in metadata}

    # Cluster label in FASTA header is the source human Ensembl gene ID.
    cluster_label_map = {
        human_gene_id: str(human_gene_id)
        for human_gene_id in ortho_clusters.keys()
        if not str(human_gene_id).startswith('_')
    }

    # Collect orthologs to fetch (deduplicate)
    seen = set()
    candidates = []
    tasks = []

    for gene_id, gene_data in ortho_clusters.items():
        if str(gene_id).startswith('_'):
            continue
        if not isinstance(gene_data, dict):
            continue

        for species, ortho_ids in gene_data.get('orthologs', {}).items():
            for o_id in ortho_ids:
                if (species, o_id) in seen:
                    continue
                seen.add((species, o_id))

                key = (o_id, species)
                if key not in metadata_lookup:
                    continue

                meta = metadata_lookup[key]
                gene_start = meta.get('start')
                gene_end = meta.get('end')

                if not gene_start or not gene_end:
                    continue

                cluster_id = cluster_label_map.get(gene_id, str(gene_id))
                candidates.append((cluster_id, o_id, species, meta))

    print(f"Computing buffers for {len(candidates)} orthologs...")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _compute_safe_buffer,
                o_id,
                species,
                meta.get('start'),
                meta.get('end'),
                deterministic,
            ): (cluster_id, o_id, species, meta)
            for cluster_id, o_id, species, meta in candidates
        }

        for future in as_completed(futures):
            cluster_id, o_id, species, meta = futures[future]
            buf_L, buf_R = future.result()
            # print(f"Calculated buffers for species {species} and gene id {o_id}: {buf_L}, {buf_R}")  # NOTE: enable for debug
            tasks.append((
                cluster_id, o_id, species, buf_L, buf_R,
                meta.get('canonical_transcript'),
                meta.get('chromosome'),
                meta.get('start'),
                meta.get('end'),
                meta.get('strand'),
            ))
                
    print(f"Fetching sequences for {len(tasks)} orthologs...")

    fasta_entries = []
    gtf_lines = []
    gff3_lines = []
    seen_ids = set()

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                _fetch_gene_data_wrapper,
                cluster_id,
                o_id,
                species,
                buf_L,
                buf_R,
                canonical_transcript_id,
                chrom,
                gene_start,
                gene_end,
                strand_int,
                return_gff3_raw=True,
            )
            for cluster_id, o_id, species, buf_L, buf_R,
                canonical_transcript_id, chrom, gene_start, gene_end, strand_int in tasks
        ]

        for future in as_completed(futures):
            result = future.result()
            if not result:
                continue

            safe_id, seq, gtf_entries, _, fasta_header, gff3_raw = result

            if safe_id in seen_ids:
                idx = 2
                base_id = safe_id
                while safe_id in seen_ids:
                    safe_id = f"{base_id}_{idx}"
                    idx += 1
            seen_ids.add(safe_id)

            is_human = 'species=homo_sapiens' in fasta_header
            fasta_entries.append((fasta_header, seq, is_human))
            gtf_lines.extend(gtf_entries)
            if gff3_raw:
                gff3_lines.extend(gff3_raw)

    fasta_dir = output_dir / "fasta"
    gtf_dir = output_dir / "gtf"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    gtf_dir.mkdir(parents=True, exist_ok=True)

    fasta_file = fasta_dir / f"{group_id}.fa"
    with open(fasta_file, 'w') as f:
        for fasta_header, seq, is_human in fasta_entries:
            wrapped_seq = _format_fasta_sequence(seq, width=80)
            if include_human or not is_human:
                f.write(f">{fasta_header}\n{wrapped_seq}\n")

    gtf_file = gtf_dir / f"{group_id}.gtf"
    with open(gtf_file, 'w') as f:
        for line in gtf_lines:
            if include_human or 'species=homo_sapiens' not in line:
                f.write(line + "\n")

    gff3_file = gtf_dir / f"{group_id}.combined.gff3"
    with open(gff3_file, 'w') as f:
        for line in gff3_lines:
            f.write(line + "\n")
    print(f"   raw combined GFF3 written to {gff3_file}")

    print(f"   {len(fasta_entries)} sequences written to {fasta_file}")
    print(f"   annotations written to {gtf_file}")

def _compute_safe_buffer(o_id: str, species: str, gene_start: int, gene_end:int, deterministic: bool = True, mu: int =5000, search_window: int =100000) -> tuple:
    """
    Calculate biologically-informed flanking buffer sizes.

    Uses an exponential distribution (mean=mu) to draw a candidate flank F,
    then caps it at the physical gap d to the nearest same-strand neighbor.
    This prevents the flanking region from overlapping neighboring genes.

    If the neighbor search fails for any reason (API error, timeout), falls
    back silently to the fixed L value so the pipeline never stalls.

    Args:
        o_id:          Ensembl stable ID of the ortholog
        species:       species name in Ensembl style (e.g. 'mus_musculus')
        gene_start:    gene start coordinate (from metadata)
        gene_end:      gene end coordinate (from metadata)
        mu:            mean of exponential distribution in bp (default: 5000)
        search_window: bp to search on each side for neighbors (default: 100000)

    """

    try:
        session = _get_thread_session()
        server = "https://rest.ensembl.org"

        # Fetch gene info (strand + chromosome) - already in metadata but
        # we need it here for the neighbor search region
        gene_resp = session.get(
            f"{server}/lookup/id/{o_id}?content-type=application/json",
            timeout=20,
        )
        if not gene_resp.ok:
            raise ValueError(f"lookup failed: {gene_resp.status_code}")
        gene_data = gene_resp.json()

        chrom  = gene_data['seq_region_name']
        strand = gene_data['strand']  # 1 or -1

        # Search for neighboring genes on the same strand within search_window
        region = f"{chrom}:{max(1, gene_start - search_window)}-{gene_end + search_window}"
        neighbors_resp = session.get(
            f"{server}/overlap/region/{species}/{region}?feature=gene"
            f";content-type=application/json",
            timeout=20,
        )
        if not neighbors_resp.ok:
            raise ValueError(f"overlap failed: {neighbors_resp.status_code}")
        neighbors = neighbors_resp.json()

        d_upstream   = float('inf')
        d_downstream = float('inf')

        for feat in neighbors:
            if feat.get('id') == o_id:
                continue
            if feat.get('strand') != strand:
                continue

            feat_start = feat['start']
            feat_end   = feat['end']

            if strand == 1:  # forward strand: upstream = left, downstream = right
                if feat_end < gene_start:
                    d_upstream   = min(d_upstream,   gene_start - feat_end - 1)
                elif feat_start > gene_end:
                    d_downstream = min(d_downstream, feat_start - gene_end - 1)
            else:             # reverse strand: upstream = right, downstream = left
                if feat_start > gene_end:
                    d_upstream   = min(d_upstream,   feat_start - gene_end - 1)
                elif feat_end < gene_start:
                    d_downstream = min(d_downstream, gene_start - feat_end - 1)

        # Deterministic exponential draw - same gene always gets same flank
        def _draw_flank(suffix) -> int:
            if deterministic:
                key  = f"{o_id}_{chrom}_{gene_start}_{suffix}"
                seed = int(hashlib.sha256(key.encode()).hexdigest(), 16) % (2 ** 32)
                rng  = np.random.default_rng(seed)
            else:
                rng  = np.random.default_rng(None)
            return rng.exponential(scale=mu)

        f_upstream   = _draw_flank("up")
        f_downstream = _draw_flank("down")

        chrom_resp = session.get(
            f"{server}/info/assembly/{species}/{chrom}?content-type=application/json",
            timeout=20,
        )
        if not chrom_resp.ok:
            raise ValueError(f"Checking chromosome length failed: {chrom_resp.status_code}")
        chrom_length = chrom_resp.json()['length']

        buf_L: int = int(max(0, min(gene_start - 1, min(d_upstream,   f_upstream))))                    
        buf_R: int = int(max(0, min(chrom_length - gene_end, min(d_downstream, f_downstream))))
                
        return buf_L, buf_R 

    except Exception as e:
        # Fallback to fixed flank so pipeline never stalls on API issues
        flank = max(0, mu)
        return flank, flank


def _fetch_gene_data_wrapper(cluster_id, o_id, species, buf_L, buf_R, canonical_transcript_id=None, chrom=None, gene_start=None, gene_end=None, strand_int=None, return_gff3_raw=False):
    """
    Fetch one ortholog: FASTA sequence with buffers and GFF3 annotation.
    Converts GFF3 to GTF, preserving original genomic coordinates.

    Notes:
    - FASTA is fetched without any masking parameter (unmasked sequence).

    Args:
        o_id: Ensembl stable ID
        species: species name
        buf_L: left buffer size
        buf_R: right buffer size
    Returns:
        tuple: (safe_id, FASTA sequence, list of GTF lines, gff_status, fasta_header, gff3_raw_lines) or None on failure
    """
    session = _get_thread_session()

    try:
        # Use coordinates from metadata - no extra lookup needed
        if not chrom or gene_start is None or gene_end is None:
            return None

        region_start = gene_start - buf_L
        region_end = gene_end + buf_R

        gff_region_start = max(1, region_start)
        gff_region_end = region_end
        strand = "+" if strand_int == 1 else "-"

        seq_url = (
            f"https://rest.ensembl.org/sequence/id/{o_id}?"
            f"expand_5prime={buf_L};expand_3prime={buf_R};content-type=text/x-fasta"
        )
        seq_resp = session.get(seq_url, timeout=30)
        if not seq_resp.ok:
            raise ValueError(f"FASTA fetch failed: HTTP {seq_resp.status_code} for {o_id}")
        seq_raw = seq_resp.text.strip().split('\n')
        if not seq_raw or not seq_raw[0].startswith('>'):
            raise ValueError(f"FASTA response malformed for {o_id}: got '{seq_raw[0][:80]}'")
        seq = ''.join(seq_raw[1:])
        if not seq:
            raise ValueError(f"FASTA sequence empty for {o_id}")

        base = f"{species}_{o_id}"
        safe_id = re.sub(r"[^A-Za-z0-9_]", "_", base)

        contig_start = max(1, region_start)
        fasta_header_extended = (
            f"cluster={cluster_id};member={o_id};species={species};"
            f"contig={chrom}:{contig_start}-{region_end};strand={strand}"
        )

        # canonical_transcript_id passed in from metadata (no extra API call needed)

        gff_url = (
            f"https://rest.ensembl.org/overlap/region/{species}/"
            f"{chrom}:{gff_region_start}..{gff_region_end}"
            f"?feature=gene;feature=transcript;feature=exon;feature=cds;content-type=text/x-gff3"
        )
        gff_resp = session.get(gff_url, timeout=30)
        if 400 <= gff_resp.status_code < 500:
            return (safe_id, seq, [], "gff_4xx", fasta_header_extended, None)
        gff_resp.raise_for_status()
        gff_text = gff_resp.text.strip()
        gff_lines = gff_text.split('\n') if gff_text else []

        gtf_entries = _gff3_to_gtf(gff_lines, fasta_header_extended, o_id, canonical_transcript_id, region_start)
        gff3_raw_lines = gff_lines if return_gff3_raw else None
        return (safe_id, seq, gtf_entries, "ok", fasta_header_extended, gff3_raw_lines)
    except Exception:
        return None


def _gff3_to_gtf(gff_lines, safe_id, gene_id, canonical_transcript_id=None, region_start=0):
    """
    Convert GFF3 feature rows to GTF format with original genomic coordinates.
    Preserves frame column from GFF3 CDS records.
    Optionally filters to only include canonical transcript features.

    Args:
        gff_lines: list of GFF3 format lines (tab-separated)
        safe_id: sanitized sequence ID for GTF output
        gene_id: Ensembl stable ID for gene_id field
        canonical_transcript_id: optional canonical transcript ID to filter by
        region_start: start coordinate of the genomic region for coordinate conversion

    Returns:
        list: GTF feature rows (tab-separated) with genomic coordinates preserved
    """
    gtf_entries = []
    gff3_feature_types = {'gene', 'exon', 'cds'}

    for line in gff_lines:
        if line.startswith('#') or not line.strip():
            continue

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            continue

        ftype = fields[2]
        ftype_lower = ftype.lower()

        if ftype_lower not in gff3_feature_types:
            continue

        source = fields[1]
        start = int(fields[3])
        end = int(fields[4])
        score = fields[5]
        strand = fields[6]
        frame = fields[7]
        attr_str = fields[8]

        attrs = parse_attributes_gff3(attr_str)

        transcript_id = attrs.get('transcript_id') or attrs.get('Parent', gene_id)
        if transcript_id.startswith('transcript:'):
            transcript_id = transcript_id.replace('transcript:', '')

        transcript_id_base = transcript_id.split('.')[0]
        canonical_id_base = canonical_transcript_id.split('.')[0] if canonical_transcript_id else None

        if canonical_id_base and transcript_id_base != canonical_id_base:
            continue

        local_start = start - region_start + 1
        local_end = end - region_start + 1

        gtf_attrs = {
            'gene_id': gene_id,
            'transcript_id': transcript_id,
        }
        if 'gene_name' in attrs:
            gtf_attrs['gene_name'] = attrs['gene_name']
        if 'biotype' in attrs:
            gtf_attrs['biotype'] = attrs['biotype']

        gtf_attr = format_attributes_gtf(gtf_attrs)
        gtf_line = f"{safe_id}\t{source}\t{ftype}\t{local_start}\t{local_end}\t{score}\t{strand}\t{frame}\t{gtf_attr}"
        gtf_entries.append(gtf_line)

    return gtf_entries