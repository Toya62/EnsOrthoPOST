#!/usr/bin/env python3
"""Species list parsing utilities for EnsOrthoPOST."""

import re
from pathlib import Path


def _to_ensembl_style(name):
    """Normalize species name to Ensembl style: lowercase_with_underscores."""
    normalized = name.strip().lower()
    normalized = re.sub(r"[^a-z0-9\s_]+", " ", normalized)
    normalized = re.sub(r"\s+", " ", normalized).strip()
    return normalized.replace(" ", "_")


def load_species_list(file_path):
    """Read one species per line, ignoring blank lines and comments."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Species file not found: {path}")

    species = []
    seen = set()
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        normalized = _to_ensembl_style(line)
        if not normalized:
            continue

        if normalized not in seen:
            seen.add(normalized)
            species.append(normalized)

    if not species:
        raise ValueError(f"No valid species found in file: {path}")

    return species

