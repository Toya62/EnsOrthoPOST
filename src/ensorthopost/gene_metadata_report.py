#!/usr/bin/env python3
"""
Fetch gene metadata (coordinates) for all orthologs
"""

import requests
import time
import json
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


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


def generate_gene_metadata(ortho_clusters, output_dir, group_id, max_workers, throttle_delay):
    """
    Fetch gene coordinates for all orthologs from Ensembl using POST batching.

    Instead of one GET per ortholog (~111k requests), we group IDs by species
    and send one POST batch per species (~139 requests). expand=1 is included
    so canonical_transcript is available for build_corpus_files.py.

    Args:
        ortho_clusters: dict of ortholog clusters
        output_dir: Path to output directory
        group_id: gene group ID
        max_workers: number of parallel workers
        throttle_delay: delay between API requests in seconds

    Returns:
        list: metadata records with coordinates and canonical_transcript
    """

    # Collect all orthologs, grouped by species for batching
    all_orthologs = []
    seen_request_pairs = set()
    species_to_ids = {}  # {species: [o_id, ...]}

    for gene_id, gene_data in ortho_clusters.items():
        if str(gene_id).startswith('_'):
            continue
        if not isinstance(gene_data, dict):
            continue

        for species, ortho_ids in gene_data.get('orthologs', {}).items():
            for o_id in ortho_ids:
                all_orthologs.append((o_id, species))
                pair = (o_id, species)
                if pair not in seen_request_pairs:
                    seen_request_pairs.add(pair)
                    species_to_ids.setdefault(species, []).append(o_id)

    unique_orthologs = list(seen_request_pairs)
    total_batches = sum(
        (len(ids) + 999) // 1000 for ids in species_to_ids.values()
    )
    print(
        f"Fetching metadata for {len(unique_orthologs)} orthologs "
        f"across {len(species_to_ids)} species "
        f"using {total_batches} POST batches (batch size ≤1000)..."
    )

    BATCH_SIZE = 1000
    POST_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
    POST_URL = "https://rest.ensembl.org/lookup/id"

    def _fetch_batch(species, ids_chunk):
        """POST one batch of up to 1000 IDs and return parsed metadata rows."""
        session = _get_thread_session()
        time.sleep(throttle_delay)
        results = []
        try:
            max_429_retries = 4
            resp = None
            for attempt in range(max_429_retries + 1):
                resp = session.post(
                    POST_URL,
                    headers=POST_HEADERS,
                    json={"ids": ids_chunk, "expand": 1},
                    timeout=60,
                )

                if resp.status_code != 429:
                    break

                if attempt == max_429_retries:
                    break

                retry_after = float(resp.headers.get("Retry-After", 0))
                backoff = min(30.0, 2.0 ** attempt)
                wait = max(retry_after, backoff)
                time.sleep(wait)

            if not resp.ok:
                # Whole batch failed - record all IDs as failed
                for o_id in ids_chunk:
                    results.append({
                        "ok": False,
                        "ensembl_id": o_id,
                        "species": species,
                        "reason": f"HTTP_{resp.status_code}_batch_error",
                    })
                return results

            data = resp.json()

            for o_id in ids_chunk:
                entry = data.get(o_id)
                if not entry:
                    results.append({
                        "ok": False,
                        "ensembl_id": o_id,
                        "species": species,
                        "reason": "null_in_batch_response",
                    })
                    continue

                start = entry.get("start")
                end = entry.get("end")
                if not start or not end:
                    results.append({
                        "ok": False,
                        "ensembl_id": o_id,
                        "species": species,
                        "reason": "missing_start_or_end",
                    })
                    continue

                # Extract canonical transcript (available because expand=1)
                canonical = entry.get("canonical_transcript")
                if not canonical:
                    for t in entry.get("Transcript", []):
                        if t.get("is_canonical") == 1:
                            canonical = t.get("id")
                            break

                results.append({
                    "ok": True,
                    "ensembl_id": o_id,
                    "species": species,
                    "start": start,
                    "end": end,
                    "strand": entry.get("strand"),
                    "chromosome": entry.get("seq_region_name"),
                    "length": end - start + 1,
                    "canonical_transcript": canonical,
                })

        except Exception as e:
            for o_id in ids_chunk:
                results.append({
                    "ok": False,
                    "ensembl_id": o_id,
                    "species": species,
                    "reason": f"request_exception:{type(e).__name__}",
                })

        return results

    # Build batch tasks: split each species into chunks of BATCH_SIZE
    batch_tasks = []
    for species, ids in species_to_ids.items():
        for i in range(0, len(ids), BATCH_SIZE):
            batch_tasks.append((species, ids[i:i + BATCH_SIZE]))

    metadata = []
    metadata_issues = []
    processed_batches = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_fetch_batch, species, chunk)
            for species, chunk in batch_tasks
        ]

        for future in as_completed(futures):
            processed_batches += 1
            if processed_batches % 10 == 0 or processed_batches == total_batches:
                print(f"  Progress: {processed_batches}/{total_batches} batches...", flush=True)

            for result in future.result():
                if result.get("ok"):
                    row = dict(result)
                    row.pop("ok", None)
                    metadata.append(row)
                else:
                    metadata_issues.append({
                        "ensembl_id": result.get("ensembl_id", ""),
                        "species": result.get("species", ""),
                        "reason": result.get("reason", "unknown_failure"),
                    })
    
    # Deduplicate metadata
    seen = set()
    unique_metadata = []
    for m in metadata:
        key = (m['ensembl_id'], m['species'])
        if key not in seen:
            seen.add(key)
            unique_metadata.append(m)

    # Deduplicate metadata issues
    seen_issues = set()
    unique_issues = []
    for issue in metadata_issues:
        key = (issue.get('ensembl_id'), issue.get('species'), issue.get('reason'))
        if key not in seen_issues:
            seen_issues.add(key)
            unique_issues.append(issue)

    total_requested = len(unique_orthologs)  # len of seen_request_pairs
    total_requested_raw = len(all_orthologs)
    duplicate_request_count = total_requested_raw - total_requested
    valid_coords = len(unique_metadata)
    failed_or_no_coords = len(unique_issues)

    # Humans are part of the ortholog set; they are included in valid_coords
    human_valid = sum(1 for m in unique_metadata if m.get('species') == 'homo_sapiens')

    print(f"   {valid_coords} genes with valid coordinates (including {human_valid} human genes)")
    print(f"   {failed_or_no_coords} orthologs without valid coordinates due to API failures or missing start/end")

    # Save to JSON file (metadata only)
    stats_dir = output_dir / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)
    metadata_file = stats_dir / f"{group_id}_ortholog_metadata.json"
    with open(metadata_file, "w", encoding="utf-8") as fh:
        json.dump(unique_metadata, fh, indent=2)
    print(f"   Saved to {metadata_file}")

    # Save metadata issues (for user-visible reasons, not only logs)
    reason_counts = {}
    for issue in unique_issues:
        reason = issue.get('reason', 'unknown_failure')
        reason_counts[reason] = reason_counts.get(reason, 0) + 1

    issues_file = stats_dir / f"{group_id}_metadata_issues.json"
    issues_report = {
        "requested_count_raw": total_requested_raw,
        "requested_count": total_requested,
        "duplicate_request_count": duplicate_request_count,
        "valid_metadata_count": valid_coords,
        "invalid_metadata_count": failed_or_no_coords,
        "reason_counts": reason_counts,
        "issues": unique_issues,
    }
    with open(issues_file, "w", encoding="utf-8") as fh:
        json.dump(issues_report, fh, indent=2)
    print(f"   Saved metadata issues to {issues_file}")

    return unique_metadata