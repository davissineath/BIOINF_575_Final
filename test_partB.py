"""Hand-worked tests for part B. Run: python test_partB.py"""

import pandas as pd

from Team7_Project import coverage_per_gene, filter_associations, person_report

# --- filter_associations ------------------------------------------------------
toy = pd.DataFrame({
    "Variant/Haplotypes": ["rs1", "rs2", "rs3", "rs4", "rs5", "CYP2D6*4", "rs6"],
    "Gene":               ["G1",  "G1",  "G1",  "G2",  "G2",  "CYP2D6",   "G3"],
    "Phenotype Category": ["Efficacy", "efficacy", "Dosage, Efficacy",
                           "Efficacy", "Toxicity", "Efficacy", "Efficacy"],
    "Significance":       ["yes", "YES", "yes", "no", "yes", "yes", "yes"],
    "Is/Is Not associated": ["Associated with", "Associated with", "Associated with",
                             "Associated with", "Associated with", "Associated with",
                             "Not associated with"],
})
kept = filter_associations(toy)
# rs4: significance no; rs5: toxicity; rs6: "Not associated with"
assert kept["Variant/Haplotypes"].tolist() == ["rs1", "rs2", "rs3", "CYP2D6*4"], kept
assert len(toy) == 7, "input must not change"
assert filter_associations(toy, "toxicity")["Variant/Haplotypes"].tolist() == ["rs5"]

# task names work too
renamed = toy.rename(columns={"Phenotype Category": "PHENOTYPE_CATEGORY",
                              "Significance": "SIGNIFICANCE",
                              "Is/Is Not associated": "ASSOCIATED"})
assert len(filter_associations(renamed)) == 4

# --- coverage_per_gene --------------------------------------------------------
pgx = pd.DataFrame({
    "Variant/Haplotypes": ["rs1", "rs1", "rs2, rs3", "rs9", "CYP2D6*4", "SLC6A4 L allele-rs25531C"],
    "Gene":               ["G1",  "G1",  "G1",       "G2",  "CYP2D6",   "SLC6A4"],
})
chip = pd.DataFrame({"dbSNP_ID": ["rs1", "rs3", "i7000", "rs100"]})
cov = coverage_per_gene(pgx, chip)
# G1: rs1 (twice, counted once), rs2, rs3 -> 3 rsIDs, rs1 and rs3 on chip
# G2: rs9 -> 1, not on chip; CYP2D6*4 and the SLC6A4 name are not rsIDs
expected = pd.DataFrame({
    "GENE_SYMBOL": ["G1", "G2"],
    "RSIDS_PHARMGKB": [3, 1],
    "RSIDS_ON_CHIP": [2, 0],
    "FRACTION_ON_CHIP": [0.667, 0.0],
})
pd.testing.assert_frame_equal(cov, expected, check_dtype=False)

# also works after A's split_column already turned the column into lists
pgx_lists = pgx.copy()
pgx_lists["Variant/Haplotypes"] = pgx_lists["Variant/Haplotypes"].str.split(",")
pd.testing.assert_frame_equal(coverage_per_gene(pgx_lists, chip), expected, check_dtype=False)

# --- person_report ------------------------------------------------------------
# stand-ins for C's functions, only for this test
def parse_alleles(text):
    return [part.strip() for part in text.split("+")]

def carrier_status(genotype, alleles):
    if genotype in ("--", ""):
        return "no call"
    if any(not set(a) <= set("ACGT") for a in alleles):
        return "unknown"
    for a in alleles:
        if (len(a) == 1 and a in genotype) or sorted(a) == sorted(genotype):
            return "yes"
    return "no"

evidence = pd.DataFrame({
    "dbSNP_ID": ["rs4244285", "rs10", "rs11", "rs12"],
    "GENE_SYMBOL": ["CYP2C19", "X", "Y", "Z"],
    "DRUG_NAME": ["clopidogrel", "a", "b", "c"],
    "LEVEL": ["1A", "1B", "2A", "2B"],
    "DIRECTION": ["decreased", "increased", "increased", "decreased"],
    "ALLELE_PharmGKB": ["AA + AG", "T", "*4", "CC"],
})
genotypes = {"rs4244285": "GA", "rs10": "CC", "rs11": "AG"}
report = person_report(evidence, genotypes, parse_alleles, carrier_status)
assert report["GENOTYPE"].tolist() == ["GA", "CC", "AG", "--"]
assert report["CARRIER"].tolist() == ["yes", "no", "unknown", "no call"]

print("all part B tests passed")
