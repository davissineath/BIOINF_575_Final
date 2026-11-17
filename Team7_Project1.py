# This file contains the functions needed for Final Project 1 in BIOINF 575
# Instructor: Cristina Mitrea Fall 2027
# Authors: Sally Alma'ani, Lauren Jeffery, Davis Sineath
###
# ---------------------------------------------------------------------------
# Part B - Sally Alma'ani
#   filter_associations   tasks 2 and 3
#   coverage_per_gene     task 3
#   plot_coverage         task 3
#   person_report         task 5 (uses C's parse_alleles and carrier_status)
#
# These functions accept either column naming:
#   - the raw PharmGKB names that read_pharmgkb / merge_chip_pharmgkb return now
#     ("Significance", "Gene", ...), or
#   - the task names ("SIGNIFICANCE", "GENE_SYMBOL", ...)
# so they keep working if the team renames columns later.
# ---------------------------------------------------------------------------

import re

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# raw PharmGKB name -> task name
PGX_TO_TASK_NAMES = {
    "Variant/Haplotypes": "dbSNP_ID",
    "Gene": "GENE_SYMBOL",
    "Drug(s)": "DRUG_NAME",
    "Phenotype Category": "PHENOTYPE_CATEGORY",
    "Significance": "SIGNIFICANCE",
    "Is/Is Not associated": "ASSOCIATED",
    "Direction of effect": "DIRECTION",
    "Notes": "NOTES",
    "Sentence": "SENTENCE",
    "Alleles": "ALLELE_PharmGKB",
}

RSID_PATTERN = re.compile(r"^rs\d+$")


def _find_column(table, *names):
    """Return the first of `names` that is a column of `table`."""
    for name in names:
        if name in table.columns:
            return name
    raise KeyError(f"none of the columns {names} is in the table")


def _clean_text(series):
    """Lower case, no surrounding spaces, missing values as ''."""
    return series.fillna("").astype(str).str.strip().str.lower()


def filter_associations(table, category="efficacy"):
    """Keep the significant, positive associations of one phenotype category.

    A row is kept when
      - Significance is "yes",
      - Is/Is Not associated is exactly "Associated with" (not "Not associated with"),
      - Phenotype Category contains `category` ("Efficacy", "Dosage, Efficacy", ...).
    Case and surrounding spaces are ignored everywhere.

    Works on the raw PharmGKB table (task 3) and on the merged table (task 2).
    Returns a new table; the input is not changed.
    """
    significance = _clean_text(table[_find_column(table, "SIGNIFICANCE", "Significance")])
    associated = _clean_text(table[_find_column(table, "ASSOCIATED", "Is/Is Not associated")])
    phenotype = _clean_text(table[_find_column(table, "PHENOTYPE_CATEGORY", "Phenotype Category")])

    keep = (
        (significance == "yes")
        & (associated == "associated with")
        & phenotype.str.contains(category.strip().lower(), regex=False)
    )
    return table.loc[keep].copy()


def _one_variant_per_row(table, column):
    """One row per variant, without surrounding spaces.

    Handles both a text column ("rs1, rs2") and a column that split_column
    already turned into lists.
    """
    out = table.copy()
    values = out[column]
    if not values.map(lambda v: isinstance(v, list)).any():
        values = values.fillna("").astype(str).str.split(",")
    out[column] = values
    out = out.explode(column)
    out[column] = out[column].fillna("").astype(str).str.strip()
    return out


def coverage_per_gene(pgx_efficacy, chip):
    """For each gene: how many of its PharmGKB rsIDs the chip measures (task 3).

    pgx_efficacy: PharmGKB annotations after filter_associations (start from
                  var_drug_ann.tsv, not from the merged table).
    chip:         the chip table, with a dbSNP_ID column.

    Haplotypes (CYP2D6*4) and other non-rsID names are dropped; each rsID is
    counted once per gene, however many publications report it. A row that lists
    several genes ("APOC1, APOE") is kept as PharmGKB lists it.

    Returns GENE_SYMBOL, RSIDS_PHARMGKB, RSIDS_ON_CHIP, FRACTION_ON_CHIP,
    largest RSIDS_PHARMGKB first (ties: more on the chip first, then by name).
    """
    variant_col = _find_column(pgx_efficacy, "Variant/Haplotypes", "dbSNP_ID")
    gene_col = _find_column(pgx_efficacy, "Gene", "GENE_SYMBOL")

    variants = _one_variant_per_row(pgx_efficacy[[gene_col, variant_col]], variant_col)
    variants[gene_col] = variants[gene_col].fillna("").astype(str).str.strip()

    is_rsid = variants[variant_col].str.match(RSID_PATTERN)
    has_gene = variants[gene_col] != ""
    pairs = (
        variants.loc[is_rsid & has_gene, [gene_col, variant_col]]
        .drop_duplicates()
        .rename(columns={gene_col: "GENE_SYMBOL", variant_col: "dbSNP_ID"})
    )

    chip_ids = set(chip["dbSNP_ID"].astype(str).str.strip())
    pairs["ON_CHIP"] = pairs["dbSNP_ID"].isin(chip_ids)

    coverage = (
        pairs.groupby("GENE_SYMBOL")
        .agg(RSIDS_PHARMGKB=("dbSNP_ID", "size"), RSIDS_ON_CHIP=("ON_CHIP", "sum"))
        .reset_index()
    )
    coverage["RSIDS_ON_CHIP"] = coverage["RSIDS_ON_CHIP"].astype(int)
    coverage["FRACTION_ON_CHIP"] = (
        coverage["RSIDS_ON_CHIP"] / coverage["RSIDS_PHARMGKB"]
    ).round(3)

    return coverage.sort_values(
        ["RSIDS_PHARMGKB", "RSIDS_ON_CHIP", "GENE_SYMBOL"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def plot_coverage(coverage, path, top_n=15):
    """Bar plot of RSIDS_PHARMGKB and RSIDS_ON_CHIP for the top_n genes (task 3).

    Horizontal bars so the gene names stay readable; the gene with the most
    rsIDs is at the top. Saves the figure to `path` and returns it.
    """
    top = coverage.sort_values("RSIDS_PHARMGKB", ascending=False).head(top_n)
    top = top.iloc[::-1]  # barh draws from the bottom up

    positions = range(len(top))
    height = 0.38
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(top) + 1.6))
    ax.barh([p + height / 2 for p in positions], top["RSIDS_PHARMGKB"],
            height=height, color="#2a78d6", label="rsIDs in PharmGKB")
    ax.barh([p - height / 2 for p in positions], top["RSIDS_ON_CHIP"],
            height=height, color="#eb6834", label="rsIDs on the 23andme chip")

    for p, (total, on_chip) in enumerate(zip(top["RSIDS_PHARMGKB"], top["RSIDS_ON_CHIP"])):
        ax.text(total + 0.3, p + height / 2, str(total), va="center", fontsize=8, color="#444444")
        ax.text(on_chip + 0.3, p - height / 2, str(on_chip), va="center", fontsize=8, color="#444444")

    ax.set_yticks(list(positions))
    ax.set_yticklabels(top["GENE_SYMBOL"])
    ax.set_xlabel("Number of rsIDs (significant efficacy associations)")
    ax.set_title(f"Chip coverage of PharmGKB efficacy variants, top {len(top)} genes")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return fig


def person_report(evidence, genotypes, parse_alleles, carrier_status):
    """The person's genotype and carrier status for each row of task 4 (task 5).

    evidence:       the task 4 table (dbSNP_ID, GENE_SYMBOL, DRUG_NAME, LEVEL,
                    DIRECTION, ALLELE_PharmGKB, ...).
    genotypes:      {rsID: genotype} from read_genotypes.
    parse_alleles,
    carrier_status: C's functions, passed in so this works as soon as they exist:
                    carrier_status(genotype, parse_alleles(allele)) -> "yes" / "no" /
                    "no call" / "unknown".

    A rsID missing from the person's file gets GENOTYPE "--" (no call).
    """
    report = evidence.copy()
    report["GENOTYPE"] = [
        str(genotypes.get(str(rsid).strip(), "--")).strip() or "--"
        for rsid in report["dbSNP_ID"]
    ]
    report["CARRIER"] = [
        carrier_status(genotype, parse_alleles(str(allele)))
        for genotype, allele in zip(report["GENOTYPE"], report["ALLELE_PharmGKB"])
    ]
    columns = ["dbSNP_ID", "GENE_SYMBOL", "DRUG_NAME", "LEVEL", "DIRECTION",
               "ALLELE_PharmGKB", "GENOTYPE", "CARRIER"]
    return report[columns]


###
#Part C: Lauren Jeffery
#read_clinical_variants
#split_column
#strongest_evidence (task 4)
#parse_alleles (task 5)
#carrier_status (task 5)
#please don't ask me to format this further

def read_clinical_variants(path):
    """
    Loading in clinical variants. Takes the path to a tsv file storing clinical data.


    Returns filtered: dataframe of columns variant (rsid), gene, evidence (level of evidence),
    and drug. 
    """
    df = pd.read_csv(
        path,
        sep="\t",
        names = ["Variant","Gene","Type","Evidence","Drug","Phenotype"],
        header = None
    )

    filtered = df[df["Variant"].str.startswith("rs", na=False)]
    filtered = split_column(filtered, "Drug")

    filtered = filtered[["Variant","Gene","Evidence","Drug"]]
    return filtered




def strongest_evidence(map_table, clinical):
    """
    Takes mapped file (output of task 2, with cols dbSNP_ID, Gene, Drug(s), Direction of effect, Sentence,
    ALLELE_PharmGKB, ALLELE_23andme) and merges it with clinical variants w/ cols Variant (rsid only),
    level of evidence, and drug.

    Outputs table with only strong evidence sorted by strength of evidence.
    Sources:
    https://www.geeksforgeeks.org/pandas/merge-join-two-dataframes-on-multiple-columns-in-pandas/
    https://www.reddit.com/r/learnpython/comments/1dle3a2/panda_filter_multiple_columns_by_multiple_values/
    https://www.tutorialspoint.com/article/python-pandas-remove-leading-and-trailing-whitespace-from-more-than-one-column

    Also used VSCode chat function to troubleshoot this function. Ghost suggestions enabled.
    """
    map_table = split_column(map_table, "Drug(s)")
    map_table["Drug(s)"] = map_table["Drug(s)"].str.strip()
    
    clinical_kept = clinical[["Variant","Evidence","Drug"]]

    merged_df = pd.merge(
        map_table,
        clinical_kept,
        left_on=["dbSNP_ID","Drug(s)"],
        right_on=["Variant","Drug"],
        how="inner"
    )
    
    filtered_merged = merged_df[merged_df["Evidence"].isin(["1A","1B","2A","2B"])]
    filtered_merged.sort_values("Evidence",ascending=False)
    
    return filtered_merged
    


def parse_alleles(input_string):
    """
    Split up alleles. I sure hope they're all formatted the same way. 
    source: https://stackoverflow.com/questions/8270092/remove-all-whitespace-in-a-string

    Ghost suggestions enabled in VSCode.
    """
    split_string = input_string.strip().replace(" ","")
    split_string = split_string.split("+")
    return split_string


def carrier_status(input_statusA, input_statusB):
    """
    Again, not super clear on this one based on the assignment, but should
    make it so the reverse allele (GA rather than AG) is still considered a 
    carrier. 

    Carrier status yields a "yes". Non-carrier status yields a "no". 
    Used VSCode chat function to troubleshoot this function. Ghost suggestions enabled. 
    """
    rev_stat_A = input_statusA[::-1]
    if rev_stat_A in input_statusB:
        return "yes"

    else:
        return "no"
