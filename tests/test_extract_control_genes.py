import ensorthopost.extract_control_genes as mod


class FakeResponse:
    def __init__(self, text="", status_code=200, headers=None, json_data=None):
        self.text = text
        self.status_code = status_code
        self.headers = headers or {}
        self._json_data = json_data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json_data


class FakeSession:
    def __init__(self, responses):
        self.responses = responses

    def mount(self, *_args, **_kwargs):
        return None

    def get(self, url, headers=None, timeout=None):
        key = url
        if key not in self.responses:
            raise AssertionError(f"Unexpected URL requested: {url}")
        result = self.responses[key]
        if isinstance(result, Exception):
            raise result
        return result


def test_extract_control_genes_download_json_success(monkeypatch):
    group_id = 141
    url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch&format=json"
    payload = [
        {
            "ensemblGeneID": "ENSG000001",
            "approvedSymbol": "GENEA",
            "locusType": "gene with protein product",
        },
        {
            "ensemblGeneID": "ENSG000002",
            "approvedSymbol": "GENEB",
            "locusType": "RNA, long non-coding",
        },
    ]
    session = FakeSession(
        {
            url: FakeResponse(
                text="[{}]",
                headers={"content-type": "application/json"},
                json_data=payload,
            )
        }
    )
    monkeypatch.setattr(mod.requests, "Session", lambda: session)

    genes = mod.extract_control_genes(group_id)

    assert genes == {"ENSG000001": "GENEA"}


def test_extract_control_genes_download_tsv_fallback(monkeypatch):
    group_id = 159
    url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch&format=json"
    tsv = (
        "Approved symbol\tEnsembl gene ID\tLocus type\n"
        "GENEA\tENSG000010\tgene with protein product\n"
        "GENEB\tENSG000011\tRNA, long non-coding\n"
    )
    session = FakeSession(
        {
            url: FakeResponse(
                text=tsv,
                headers={"content-type": "text/plain"},
                json_data=None,
            )
        }
    )
    monkeypatch.setattr(mod.requests, "Session", lambda: session)

    genes = mod.extract_control_genes(group_id)

    assert genes == {"ENSG000010": "GENEA"}


def test_extract_control_genes_rest_api_fallback(monkeypatch):
    group_id = 595
    download_url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch&format=json"
    hgnc_url = f"https://rest.genenames.org/search/gene_group_id:{group_id}"
    ensa_url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/GENEA?expand=0"
    ensb_url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/GENEB?expand=0"

    session = FakeSession(
        {
            download_url: RuntimeError("download endpoint unavailable"),
            hgnc_url: FakeResponse(
                status_code=200,
                json_data={"response": {"docs": [{"symbol": "GENEA"}, {"symbol": "GENEB"}]}}
            ),
            ensa_url: FakeResponse(status_code=200, json_data={"id": "ENSG000020", "biotype": "protein_coding"}),
            ensb_url: FakeResponse(status_code=200, json_data={"id": "ENSG000021", "biotype": "lncRNA"}),
        }
    )
    monkeypatch.setattr(mod.requests, "Session", lambda: session)

    genes = mod.extract_control_genes(group_id)

    assert genes == {"ENSG000020": "GENEA"}


def test_extract_control_genes_returns_empty_when_all_paths_fail(monkeypatch):
    group_id = 999999
    download_url = f"https://www.genenames.org/cgi-bin/genegroup/download?id={group_id}&type=branch&format=json"
    hgnc_url = f"https://rest.genenames.org/search/gene_group_id:{group_id}"

    session = FakeSession(
        {
            download_url: RuntimeError("download failed"),
            hgnc_url: FakeResponse(status_code=500, json_data={}),
        }
    )
    monkeypatch.setattr(mod.requests, "Session", lambda: session)

    genes = mod.extract_control_genes(group_id)

    assert genes == {}
