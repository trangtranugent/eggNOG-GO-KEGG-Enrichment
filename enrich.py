#!/usr/bin/env python3
"""Reusable eggNOG-mapper GO/KEGG ORA database builder and DEG enrichment CLI."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

MISSING = {"", "-", "NA", "N/A", "nan", "None", "null"}
ROOT_GO = {"GO:0008150", "GO:0003674", "GO:0005575"}
GO_RE = re.compile(r"^GO:\d{7}$")
KEGG_RE = re.compile(r"^(?:ko|map)\d{5}$")


def clean_str(v):
    if pd.isna(v):
        return ""
    value = str(v).strip()
    return "" if value in MISSING else value


def split_terms(raw):
    return {t.strip() for t in clean_str(raw).split(",") if t.strip() and t.strip() not in MISSING}


def normalize_gene(value, rule):
    """Return the original ID, or remove a user-chosen terminal isoform pattern.

    A pattern is selected against actual RNA-seq IDs, not assumed from species.
    """
    gene = clean_str(value)
    if rule is True:  # Legacy support for old databases
        rule = r"_T\d+$"
    if isinstance(rule, str) and rule:
        gene = re.sub(rule, "", gene, flags=re.IGNORECASE)
    return gene


def read_emapper(path):
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Annotation file not found: {path}")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        # eggNOG exports sometimes include metadata rows above the actual header.
        preview = pd.read_excel(path, header=None, dtype=str, nrows=150)
        matches = [i for i in range(len(preview))
                   if str(preview.iloc[i, 0]).strip().lstrip("#").lower() == "query"]
        if not matches:
            raise ValueError("Excel eggNOG file must contain a 'query' header (first column).")
        df = pd.read_excel(path, header=matches[0], dtype=str)
    else:
        # eggNOG-mapper v2 typically includes many '##' comment lines and a '#query' header.
        # Some exported tables instead use a plain 'query' header.
        with path.open("r", encoding="utf-8-sig") as handle:
            first_header = None
            for idx, line in enumerate(handle):
                trimmed = line.lstrip().lstrip("#").strip()
                first = trimmed.split("\t", 1)[0].strip().lower()
                if first == "query":
                    first_header = idx
                    break
        if first_header is None:
            raise ValueError("Could not locate a 'query' or '#query' tab-delimited header in eggNOG file")
        df = pd.read_csv(path, sep="\t", skiprows=first_header, dtype=str, low_memory=False,
                         comment=None, on_bad_lines="warn")
    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]
    aliases = {c.casefold(): c for c in df.columns}
    for target in ("query", "GOs", "KEGG_Pathway", "KEGG_ko"):
        if target not in df.columns and target.casefold() in aliases:
            df = df.rename(columns={aliases[target.casefold()]: target})
    if "query" not in df:
        raise ValueError("No query column found; check the eggNOG annotation table")
    for name in ("GOs", "KEGG_Pathway", "KEGG_ko"):
        if name not in df:
            df[name] = ""
    return df[df["query"].map(clean_str).ne("")].copy()


def read_go_obo(path):
    terms = {}
    if not path:
        return terms
    current = {}
    def flush():
        if current.get("id") and current.get("is_obsolete") != "true":
            id_ = current["id"]
            terms[id_] = (current.get("name", ""), current.get("namespace", ""))
            for alias in current.get("alt_id", []):
                terms[alias] = (current.get("name", ""), current.get("namespace", ""))
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                flush()
                current = {"alt_id": []}
            elif line.startswith("["):
                flush()
                current = {}
            elif ": " in line and current:
                key, value = line.split(": ", 1)
                if key == "alt_id":
                    current["alt_id"].append(value)
                elif key in ("id", "name", "namespace", "is_obsolete"):
                    current[key] = value
    flush()
    return terms


def read_kegg_names(path):
    names = {}
    if path:
        with open(path, encoding="utf-8-sig", newline="") as f:
            for row in csv.reader(f, delimiter="\t"):
                if len(row) >= 2:
                    key = row[0].strip()
                    if key.lower() in {"term", "pathway_id", "id"}:
                        continue
                    key = key.removeprefix("path:")
                    if re.fullmatch(r"(?:map|ko)\d{5}", key):
                        names["map" + key[-5:]] = row[1].strip()
    return names


def build(args):
    df = read_emapper(args.emapper)
    go_names = read_go_obo(args.go_obo)
    kegg_names = read_kegg_names(args.kegg_names)
    dbpath = Path(args.outdir) / "annotation.sqlite"
    dbpath.parent.mkdir(parents=True, exist_ok=True)
    gene_terms = set()
    gene_kos = set()
    genes = set()
    malformed_go, malformed_pathways = set(), set()
    for record in df.itertuples(index=False):
        values = record._asdict()
        gene = normalize_gene(values["query"], getattr(args, "id_regex", None) or (r"_T\d+$" if getattr(args, "legacy_collapse", False) else None))
        genes.add(gene)
        for go in split_terms(values["GOs"]):
            if GO_RE.fullmatch(go):
                if not args.include_go_roots and go in ROOT_GO:
                    continue
                gene_terms.add((gene, "GO", go))
            else:
                malformed_go.add(go)
        for term in split_terms(values["KEGG_Pathway"]):
            term = term.removeprefix("path:")
            if KEGG_RE.fullmatch(term):
                gene_terms.add((gene, "KEGG", "map" + term[-5:]))
            else:
                malformed_pathways.add(term)
        for ko in split_terms(values["KEGG_ko"]):
            ko = ko.removeprefix("ko:")
            if re.fullmatch(r"K\d{5}", ko):
                gene_kos.add((gene, ko))
    if dbpath.exists() and not args.force:
        raise ValueError(f"Database exists: {dbpath}. Pass --force to replace it.")
    if dbpath.exists():
        dbpath.unlink()
    with sqlite3.connect(dbpath) as con:
        con.executescript("""
            CREATE TABLE genes (gene_id TEXT PRIMARY KEY);
            CREATE TABLE terms (category TEXT NOT NULL, term_id TEXT NOT NULL,
              term_name TEXT NOT NULL, namespace TEXT NOT NULL, PRIMARY KEY (category, term_id));
            CREATE TABLE gene_terms (gene_id TEXT NOT NULL, category TEXT NOT NULL, term_id TEXT NOT NULL,
              PRIMARY KEY (gene_id, category, term_id));
            CREATE INDEX idx_gene_terms_category_term ON gene_terms(category, term_id);
            CREATE INDEX idx_gene_terms_gene ON gene_terms(gene_id);
            CREATE TABLE gene_kos (gene_id TEXT NOT NULL, ko_id TEXT NOT NULL, PRIMARY KEY(gene_id, ko_id));
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        """)
        con.executemany("INSERT INTO genes VALUES (?)", [(g,) for g in sorted(genes)])
        term_ids = {(cat, term) for _, cat, term in gene_terms}
        terms = []
        for cat, term in sorted(term_ids):
            name, namespace = go_names.get(term, ("", "")) if cat == "GO" else (kegg_names.get(term, ""), "")
            terms.append((cat, term, name, namespace))
        con.executemany("INSERT INTO terms VALUES (?,?,?,?)", terms)
        con.executemany("INSERT INTO gene_terms VALUES (?,?,?)", sorted(gene_terms))
        con.executemany("INSERT INTO gene_kos VALUES (?,?)", sorted(gene_kos))
        metadata = dict(source_emapper=str(Path(args.emapper).resolve()),
                        date_utc=datetime.now(timezone.utc).isoformat(),
                        gene_id_policy="regex:" + args.id_regex if getattr(args, "id_regex", None) else "keep",
                        kegg_pathway_policy="normalize koXXXXX/mapXXXXX -> mapXXXXX",
                        go_roots_included=str(args.include_go_roots))
        con.executemany("INSERT INTO metadata VALUES (?,?)", list(metadata.items()))
    output = Path(args.outdir)
    pd.DataFrame(sorted(genes), columns=["gene_id"]).to_csv(output / "all_genes.tsv", sep="\t", index=False)
    for category, filename in (("GO", "gene2go.tsv"), ("KEGG", "gene2kegg.tsv")):
        rows = [(g, t) for g, c, t in sorted(gene_terms) if c == category]
        pd.DataFrame(rows, columns=["gene_id", "term_id"]).to_csv(output / filename, sep="\t", index=False)
    pd.DataFrame(terms, columns=["category", "term_id", "term_name", "namespace"]).to_csv(
        output / "term_metadata.tsv", sep="\t", index=False)
    pd.DataFrame(sorted(gene_kos), columns=["gene_id", "ko_id"]).to_csv(
        output / "gene2ko.tsv", sep="\t", index=False)
    stats = dict(input_rows=len(df), unique_genes=len(genes), go_gene_term_links=sum(c == "GO" for _, c, _ in gene_terms),
                 kegg_gene_term_links=sum(c == "KEGG" for _, c, _ in gene_terms),
                 go_terms=sum(c == "GO" for c, _ in term_ids), kegg_pathways=sum(c == "KEGG" for c, _ in term_ids),
                 ko_links=len(gene_kos), invalid_go_tokens=sorted(malformed_go)[:15],
                 invalid_kegg_tokens=sorted(malformed_pathways)[:15])
    (output / "build_summary.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(f"Database: {dbpath}\n" + json.dumps(stats, indent=2))


def read_gene_table(path, requested_column, collapse, p_col=None, p_max=None, fc_col=None, direction="all"):
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, dtype=str)
    else:
        # Detect CSV/TSV; a single-column gene list is also a valid input.
        with path.open("r", encoding="utf-8-sig") as handle:
            first = next((line for line in handle if line.strip() and not line.startswith("#")), "")
        sep = "\t" if "\t" in first else "," if "," in first else "\t"
        df = pd.read_csv(path, sep=sep, dtype=str, comment="#")
        if len(df.columns) == 1 and not requested_column and not p_col and not fc_col:
            # A one-column list without a header: accept both with/without 'gene_id' header.
            name = str(df.columns[0])
            if name.casefold() not in {"gene", "gene_id", "geneid", "query", "id", "feature", "feature_id"}:
                df = pd.read_csv(path, sep="\t", header=None, dtype=str, names=["gene_id"], comment="#")
    if requested_column:
        if requested_column not in df:
            raise ValueError(f"Column '{requested_column}' not present in {path}: {list(df.columns)}")
        gene_col = requested_column
    else:
        choices = ("gene_id", "GeneID", "gene", "Gene", "query", "Geneid", "feature_id", "Feature", "id", "ID")
        gene_col = next((c for c in choices if c in df.columns), None)
        if not gene_col:
            if len(df.columns) == 1:
                gene_col = df.columns[0]
            else:
                raise ValueError(f"Specify --gene-column. Columns in {path}: {list(df.columns)}")
    if p_col:
        if p_col not in df:
            raise ValueError(f"Missing p-adjusted column: {p_col}")
        df = df[pd.to_numeric(df[p_col], errors="coerce") <= p_max]
    if direction != "all":
        if not fc_col or fc_col not in df:
            raise ValueError("--direction up/down requires --lfc-column matching a column in the DEG table")
        fcs = pd.to_numeric(df[fc_col], errors="coerce")
        df = df[fcs.gt(0) if direction == "up" else fcs.lt(0)]
    return {g for g in (normalize_gene(v, collapse) for v in df[gene_col]) if g}


def bh_adjust(pvalues):
    # Benjamini-Hochberg correction, independently for GO and KEGG tests
    n = len(pvalues)
    if not n:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    result = [1.0] * n
    q = 1.0
    for rank_index in range(n - 1, -1, -1):
        i = order[rank_index]
        q = min(q, pvalues[i] * n / (rank_index + 1))
        result[i] = q
    return result


def enrich(args):
    dbpath = Path(args.database)
    if not dbpath.is_file():
        raise ValueError(f"Database not found: {dbpath}")
    with sqlite3.connect(dbpath) as con:
        genes_db = {r[0] for r in con.execute("SELECT gene_id FROM genes")}
        gene_policy = dict(con.execute("SELECT key,value FROM metadata"))["gene_id_policy"]
        collapse = gene_policy.removeprefix("regex:") if gene_policy.startswith("regex:") else gene_policy.startswith("strip")
        annotations = defaultdict(lambda: defaultdict(set))
        for gene, cat, term in con.execute("SELECT gene_id,category,term_id FROM gene_terms"):
            annotations[cat][term].add(gene)
        term_info = {(cat, term): (name, namespace) for cat, term, name, namespace in
                     con.execute("SELECT category,term_id,term_name,namespace FROM terms")}
    de_genes = read_gene_table(args.deg, args.gene_column, collapse,
                               args.padj_column, args.padj_cutoff, args.lfc_column, args.direction)
    if args.background:
        background = read_gene_table(args.background, args.background_column, collapse)
    else:
        background = genes_db.copy()
        print("WARNING: using all genes from eggNOG annotation as background. For RNA-seq, supply --background with all genes tested for differential expression.", file=sys.stderr)
    outside_deg = len(de_genes - background)
    universe = background & genes_db if args.universe == "annotated" else background
    targets = de_genes & universe
    if not universe:
        raise ValueError("No background genes matched the annotation database. Check gene IDs/transcript handling.")
    if not targets:
        raise ValueError("No DEG genes matched the selected universe. Check --gene-column and gene IDs.")
    if args.universe == "annotated" and (background - genes_db):
        print(f"NOTE: {len(background - genes_db)} background genes have no eggNOG row and were excluded (annotated universe).", file=sys.stderr)
    results = {}
    for cat in ("GO", "KEGG"):
        # For the annotated mode, constrain to genes with at least one annotation in this category.
        cat_genes = set().union(*annotations[cat].values()) if annotations[cat] else set()
        pop = universe & cat_genes if args.universe == "annotated" else universe
        sample = targets & pop
        M, N = len(pop), len(sample)
        rows = []
        if M and N:
            for term, annotated in annotations[cat].items():
                members = pop & annotated
                K = len(members)
                if not args.min_size <= K <= args.max_size:
                    continue
                overlap = sample & members
                k = len(overlap)
                if k < args.min_overlap:
                    continue
                p = float(hypergeom.sf(k - 1, M, K, N))
                name, namespace = term_info.get((cat, term), ("", ""))
                rows.append(dict(category=cat, term_id=term, description=name, namespace=namespace,
                                 background_size=M, input_size=N, term_size=K, overlap_count=k,
                                 enrichment_ratio=(k / N) / (K / M), pvalue=p,
                                 overlapping_genes=",".join(sorted(overlap))))
        # Correct over every term passing background-size filters, not just nonzero hits.
        # The number of hypotheses is critical to FDR validity.
        if M and N:
            number_tests = sum(args.min_size <= len(pop & gs) <= args.max_size for gs in annotations[cat].values())
        else:
            number_tests = 0
        if rows:
            # Equivalent to BH including tests with zero overlap (p=1), because all retained tests have p<1.
            ps = [r["pvalue"] for r in rows]
            m = number_tests
            order = sorted(range(len(ps)), key=lambda i: ps[i])
            q = 1.0
            adjusted = [1.0] * len(ps)
            for rank_idx in range(len(order) - 1, -1, -1):
                idx = order[rank_idx]
                q = min(q, ps[idx] * m / (rank_idx + 1))
                adjusted[idx] = q
            for row, qval in zip(rows, adjusted):
                row["padj"] = qval
        columns = ["category", "term_id", "description", "namespace", "background_size", "input_size",
                   "term_size", "overlap_count", "enrichment_ratio", "pvalue", "padj", "overlapping_genes"]
        frame = pd.DataFrame(rows, columns=columns)
        if len(frame):
            frame = frame.sort_values(["padj", "pvalue", "overlap_count"], ascending=[True, True, False])
        results[cat] = dict(df=frame, tested_terms=number_tests, input_size=N, background_size=M)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for cat, info in results.items():
        frame = info["df"]
        frame.to_csv(outdir / f"{cat.lower()}_all_results.tsv", sep="\t", index=False)
        frame[frame.padj <= args.fdr].to_csv(outdir / f"{cat.lower()}_significant.tsv", sep="\t", index=False)
    summary = dict(deg_file=str(args.deg), background_file=str(args.background) if args.background else "ALL_ANNOTATED_GENES_FALLBACK",
                   universe_policy=args.universe, input_unique_degs=len(de_genes),
                   background_unique_genes=len(background), degs_outside_background=outside_deg,
                   degs_matching_selected_universe=len(targets),
                   categories={cat: {"background_size": info["background_size"], "deg_size": info["input_size"],
                                     "terms_tested": info["tested_terms"],
                                     "terms_with_reported_overlap": len(info["df"]),
                                     "significant_terms": int((info["df"].padj <= args.fdr).sum())}
                               for cat, info in results.items()})
    (outdir / "enrichment_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Results written to {outdir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    build_parser = sub.add_parser("build", help="Build a reusable SQLite and TSV annotation database")
    build_parser.add_argument("--emapper", required=True, help="eggNOG-mapper .annotations/.tsv/.xlsx file")
    build_parser.add_argument("--outdir", required=True)
    build_parser.add_argument("--go-obo", help="Optional go-basic.obo for GO names and namespaces")
    build_parser.add_argument("--kegg-names", help="Optional tab-separated pathway ID and description file")
    build_parser.add_argument("--id-regex", help="Optional terminal isoform suffix regex to remove from annotation IDs (default: preserve IDs)")
    build_parser.add_argument("--include-go-roots", action="store_true", help="Retain GO ontology root terms")
    build_parser.add_argument("--force", action="store_true", help="Replace existing database")
    build_parser.set_defaults(func=build)
    run = sub.add_parser("enrich", help="GO and KEGG pathway over-representation analysis")
    run.add_argument("--database", required=True, help="Path to annotation.sqlite")
    run.add_argument("--deg", required=True, help="CSV, TSV or one-per-line gene list")
    run.add_argument("--background", help="CSV, TSV or gene list of ALL genes tested in RNA-seq (strongly recommended)")
    run.add_argument("--gene-column", help="Gene identifier column in DEG table")
    run.add_argument("--background-column", help="Gene identifier column in background table")
    run.add_argument("--padj-column", help="Optional DEG padj column; without this, uses all rows")
    run.add_argument("--padj-cutoff", type=float, default=0.05)
    run.add_argument("--lfc-column", help="Column with log2 fold change")
    run.add_argument("--direction", choices=["all", "up", "down"], default="all")
    run.add_argument("--universe", choices=["annotated", "all"], default="annotated",
                     help="Annotated universe filters to genes with GO/KEGG annotation, separately per category")
    run.add_argument("--min-size", type=int, default=3)
    run.add_argument("--max-size", type=int, default=2000)
    run.add_argument("--min-overlap", type=int, default=1)
    run.add_argument("--fdr", type=float, default=0.05)
    run.add_argument("--outdir", required=True)
    run.set_defaults(func=enrich)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
