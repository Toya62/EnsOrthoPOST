from pathlib import Path

import pytest

from ensorthopost.helper_define_species import load_species_list


def test_load_species_list_normalizes_and_deduplicates(tmp_path: Path):
    species_file = tmp_path / "species.txt"
    species_file.write_text(
        """
# Comment line should be ignored
Homo sapiens
homo_sapiens
Canis lupus familiaris
Canis lupus familiaris
Pan troglodytes!

""",
        encoding="utf-8",
    )

    species = load_species_list(species_file)

    assert species == [
        "homo_sapiens",
        "canis_lupus_familiaris",
        "pan_troglodytes",
    ]


def test_load_species_list_raises_for_missing_file(tmp_path: Path):
    missing = tmp_path / "does_not_exist.txt"

    with pytest.raises(FileNotFoundError):
        load_species_list(missing)


def test_load_species_list_raises_for_empty_or_comment_only_file(tmp_path: Path):
    species_file = tmp_path / "species.txt"
    species_file.write_text("\n# only comments\n   \n", encoding="utf-8")

    with pytest.raises(ValueError):
        load_species_list(species_file)
