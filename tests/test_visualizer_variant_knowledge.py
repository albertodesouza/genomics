"""Tests for the known-variant annotation of the gene-products page (genomics.visualizer.variant_knowledge).

Ensembl's answers are canned: VEP's co-located variants and the variation endpoint's phenotypes.
"""
import json

import pytest

from genomics.visualizer.remote import RemoteError
from genomics.visualizer.variant_knowledge import VariantKnowledge, _known_record, _pheno_name, variant_key

VEP_ITEM = {"input": "16 89919709 . C T . . .", "colocated_variants": [
    {"id": "CM981238", "allele_string": "HGMD_MUTATION", "phenotype_or_disease": 1},
    {"id": "rs1805007", "allele_string": "C/T", "phenotype_or_disease": 1, "clin_sig": ["benign", "risk_factor", "association"],
     "pubmed": [9571181, 17952075], "var_synonyms": {"ClinVar": ["RCV000015387", "VCV000014312"], "OMIM": [155555.0004], "UniProt": ["VAR_008522"]},
     "frequencies": {"T": {"gnomadg": 0.046, "eur": 0.072, "af": 0.02}}},
]}
PHENOTYPES = {"name": "rs1805007", "clinical_significance": ["benign"], "phenotypes": [
    {"trait": "SKIN/HAIR/EYE PIGMENTATION 2, RED HAIR/FAIR SKIN", "source": "ClinVar", "risk_allele": "T"},
    {"trait": "Increased analgesia from kappa-opioid receptor agonist, female-specific", "source": "ClinVar"},
    {"trait": "ClinVar: phenotype not specified", "source": "ClinVar"},
    {"trait": "Red vs. brown/black hair color", "source": "NHGRI-EBI GWAS catalog", "pvalue": "1e-300", "study": "PMID:30531825"},
    {"trait": "Basal cell carcinoma PheCode 172.21", "source": "NHGRI-EBI GWAS catalog", "pvalue": "5e-123", "study": "PMID:39024449"},
    {"trait": "Basal cell carcinoma", "source": "NHGRI-EBI GWAS catalog", "pvalue": "4e-17", "study": "PMID:21700618"},
]}


class FakeRemote:
    offline = False

    def __init__(self, answers):
        self.answers = answers

    def get_json(self, url, ttl=None, params=None):
        for key, value in self.answers.items():
            if url.endswith(key):
                return value
        raise RemoteError("not found")


def test_known_record_picks_the_dbsnp_variant_of_the_same_allele():
    rec = _known_record(VEP_ITEM, "C", "T")
    assert rec["rsid"] == "rs1805007" and rec["same_allele"] and rec["known"]
    assert rec["clinvar"] == ["VCV000014312"] and rec["omim"] == ["155555.0004"] and rec["uniprot"] == ["VAR_008522"]
    assert rec["pubmed_count"] == 2 and rec["frequency"]["gnomad_genomes"] == 0.046 and rec["clinical_significance"][1] == "risk factor"
    assert _known_record({"colocated_variants": [{"id": "COSV1", "allele_string": "COSMIC_MUTATION"}]}, "C", "T") is None
    assert _known_record(None, "C", "T") is None
    assert _pheno_name("Basal cell carcinoma PheCode 172.21") == "Basal cell carcinoma"


def test_annotate_batches_vep_and_reads_traits(tmp_path, monkeypatch):
    knowledge = VariantKnowledge(FakeRemote({"rs1805007": PHENOTYPES}), tmp_path)
    sent = []

    def fake_post(url, payload):
        sent.append(payload)
        return [VEP_ITEM, {"input": "16 89919000 . A G . . .", "colocated_variants": []}]

    monkeypatch.setattr(knowledge, "_post", fake_post)
    out = knowledge.annotate("chr16", [(89919709, "C", "T"), (89919000, "A", "G"), (89919001, "A", "<DEL>")])
    assert out["available"] and set(out["variants"]) == {variant_key(89919709, "C", "T")}
    assert sent[0]["variants"] == ["16 89919709 . C T . . .", "16 89919000 . A G . . ."] and sent[0]["pubmed"] == 1  # symbolic alleles are not sent
    rec = out["variants"]["89919709:C:T"]
    traits = [t["trait"] for t in rec["traits"]]
    # ClinVar conditions come first, then GWAS traits by p-value
    assert set(traits[:2]) == {"SKIN/HAIR/EYE PIGMENTATION 2, RED HAIR/FAIR SKIN", "Increased analgesia from kappa-opioid receptor agonist, female-specific"}
    assert traits[2] == "Red vs. brown/black hair color"
    assert "ClinVar: phenotype not specified" not in traits and traits.count("Basal cell carcinoma") == 1  # PheCode variants merged
    bcc = next(t for t in rec["traits"] if t["trait"] == "Basal cell carcinoma")
    assert bcc["best_p"] == pytest.approx(5e-123) and len(bcc["studies"]) == 2


def test_offline_is_reported_and_responses_are_cached(tmp_path):
    offline = FakeRemote({})
    offline.offline = True
    knowledge = VariantKnowledge(offline, tmp_path)
    out = knowledge.annotate("16", [(1, "A", "G")])
    assert out["available"] is False and "no-remote" in out["error"]
    # a cached VEP answer is served even offline
    import hashlib

    from genomics.visualizer.variant_knowledge import VEP_REGION

    body = json.dumps({"variants": ["16 89919709 . C T . . ."], "pubmed": 1, "var_synonyms": 1}, sort_keys=True).encode()
    (tmp_path / "variants").mkdir()
    (tmp_path / "variants" / f"{hashlib.sha1(VEP_REGION.encode() + body).hexdigest()}.json").write_text(json.dumps([VEP_ITEM]), encoding="utf-8")
    out = knowledge.annotate("16", [(89919709, "C", "T")], with_phenotypes=False)
    assert out["available"] and out["variants"]["89919709:C:T"]["rsid"] == "rs1805007"
