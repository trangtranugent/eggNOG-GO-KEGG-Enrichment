#!/usr/bin/env python3
"""One-command workflow: add an eggNOG annotation and *.up.txt/*.down.txt DEG lists."""
from __future__ import annotations
import argparse
import csv
import json
import re
import shutil
import urllib.request
import urllib.error
from datetime import datetime, timezone
import sys
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
from enrich import build, enrich, read_emapper, read_gene_table, normalize_gene

ROOT = Path(__file__).resolve().parent

# Conservative, terminal isoform patterns. Others are supported via --id-regex.
ISOFORM_PATTERNS = [
    ("underscore_T", r"_T\d+$"),       # Lm7210a001g000030_T001
    ("dot_t", r"\.t\d+$"),             # gene123.t1
    ("dot_T", r"\.T\d+$"),             # gene123.T1
    ("hyphen_T", r"-T\d+$"),          # gene123-T1
    ("dot_p", r"\.p\d+$"),             # gene123.p1
    ("underscore_i", r"_i\d+$"),      # Trinity isoform naming
]


def choose_id_rule(annot_ids, input_ids, mode, custom_regex=None):
    """Choose one consistent ID rule for annotation, DEG and background IDs."""
    if mode == "exact":
        return None, {"policy": "exact", "exact_matches": len(annot_ids & input_ids)}
    if mode == "regex":
        if not custom_regex:
            raise ValueError("--id-mode regex requires --id-regex")
        import re
        re.compile(custom_regex)
        return custom_regex, {"policy": "custom regex", "pattern": custom_regex}
    exact = len(annot_ids & input_ids)
    candidates = [(exact, "exact", None)]
    for name, pat in ISOFORM_PATTERNS:
        normalized_annot = {normalize_gene(g, pat) for g in annot_ids}
        normalized_input = {normalize_gene(g, pat) for g in input_ids}
        candidates.append((len(normalized_annot & normalized_input), name, pat))
    # Only normalize when it actually increases the number of matched identifiers.
    count, name, pat = max(candidates, key=lambda item: item[0])
    return pat, {"policy": name, "exact_matches": exact, "matched_after_rule": count,
                 "pattern": pat, "all_candidate_matches": {n: c for c, n, _ in candidates}}


def find_annotation(input_dir: Path) -> Path:
    names = ["emapper.annotations", "emapper.annotations.tsv", "emapper.annotations.txt", "emapper.annotations.xlsx", "emapper.annotations.xls", "emapper.annotations.csv"]
    files = [input_dir / n for n in names if (input_dir / n).is_file()]
    if len(files) != 1:
        raise ValueError(f"Expected exactly ONE annotation file in {input_dir} named emapper.annotations[.tsv/.txt/.xlsx/.xls/.csv]. Found: {[p.name for p in files]}")
    return files[0]

def discover_deg_pairs(input_dir: Path):
    """Identify any <comparison>.up.txt/<comparison>.down.txt files, case-insensitively.

    Files named DEGs.up.txt and DEGs.down.txt are treated as the unnamed/default
    comparison for backward compatibility. An all list is ALWAYS computed as the
    unique union of the up and down genes; supplied *.all.txt files are ignored.
    """
    pairs = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        m = re.fullmatch(r"(.+)\.(up|down)\.txt", path.name, re.IGNORECASE)
        if not m:
            continue
        stem, direction = m.group(1), m.group(2).lower()
        if not stem.strip():
            continue
        comparison = 'default' if stem.lower() == 'degs' else (stem[5:] if stem.lower().startswith('degs.') else stem)
        # Avoid writing outside output (comparison becomes a path component).
        if comparison in ('', '.', '..') or '/' in comparison or '\\' in comparison:
            raise ValueError(f'Invalid comparison name derived from {path.name}: {comparison}')
        key=(comparison, direction)
        if key in pairs:
            raise ValueError(f'Duplicate DEG input for {comparison}/{direction}: {pairs[key]} and {path}')
        pairs[key]=path
    return pairs


def write_gene_list(path: Path, genes):
    """Write an unambiguous, one-column gene-ID TSV usable by enrich.py."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('gene_id\n' + ''.join(g + '\n' for g in sorted(genes)), encoding='utf-8')


RESOURCE_URLS = {
    "go": "https://purl.obolibrary.org/obo/go/go-basic.obo",
    "kegg": "https://rest.kegg.jp/list/pathway",
}


def get_resource(path: Path, url: str, validate, label: str, offline=False):
    """Use bundled/cached reference if valid, otherwise download once and cache."""
    if path.is_file() and validate(path):
        print(f"{label}: using cached {path}", flush=True)
        return path
    if offline:
        raise RuntimeError(f"{label} reference missing/invalid: {path}. "
                           "Run without --offline once to download it, or place the reference there.")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.download')
    print(f"{label}: fetching {url}", flush=True)
    try:
        req=urllib.request.Request(url, headers={"User-Agent":"TrangTran-eggnog-enrichment/1.1 (academic research)"})
        with urllib.request.urlopen(req,timeout=90) as response, temp.open('wb') as fout:
            shutil.copyfileobj(response, fout)
        if not validate(temp):
            raise ValueError('download did not pass basic format checks')
        temp.replace(path)
        return path
    except Exception as exc:
        temp.unlink(missing_ok=True)
        raise RuntimeError(f"Cannot obtain the {label} reference from {url}: {exc}. "
                           f"Download it manually to {path} and rerun, or try again with Internet access.") from exc


def valid_go(path):
    # A full ontology is ~tens of MB; disallow accidentally using a demo subset.
    if path.stat().st_size < 1_000_000:
        return False
    with path.open('r',encoding='utf-8') as f:
        head=''.join(f.readline() for _ in range(120))
    return 'format-version:' in head and 'ontology: go' in head


def valid_kegg(path):
    # KEGG REST list/pathway: path:map00010<TAB>Glycolysis / Gluconeogenesis
    rows=path.read_text(encoding='utf-8').splitlines()
    return len(rows)>100 and sum(bool(re.match(r"^(?:path:)?(?:map|ko)\d{5}\t\S", line)) for line in rows)>100


def resource_manifest(output: Path, go: Path, kegg: Path):
    import hashlib
    def info(path,source):
        return {'file': str(path), 'source':source,
                'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'bytes':path.stat().st_size}
    (output/'resource_manifest.json').write_text(json.dumps({
        'go':info(go,RESOURCE_URLS['go']), 'kegg':info(kegg,RESOURCE_URLS['kegg']),
        'recorded_at_utc':datetime.now(timezone.utc).isoformat()
    },indent=2),encoding='utf-8')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', default=str(ROOT/'input'), help='Directory with emapper.annotations and *.up.txt / *.down.txt DEG tables')
    ap.add_argument('--output', default=str(ROOT/'output'), help='Destination for database and enrichment tables')
    ap.add_argument('--background', help='Optional list of all RNA-seq genes eligible for DEG testing; if absent looks for input/background.txt')
    ap.add_argument('--go-obo', help='Optional GO OBO file to populate term names and GO categories')
    ap.add_argument('--kegg-names', help='Optional KEGG ID-to-name TSV')
    ap.add_argument('--offline',action='store_true',help='Never download resources; require already cached GO and KEGG metadata')
    ap.add_argument('--id-mode', choices=['auto','exact','regex'], default='auto',
                    help='auto matches IDs and only removes recognized isoform endings if this improves coverage')
    ap.add_argument('--id-regex', help=r'Custom terminal isoform pattern, e.g. _isoform\d+$ (with --id-mode regex)')
    ap.add_argument('--universe', choices=['annotated','all'], default='annotated')
    ap.add_argument('--min-size', type=int, default=3)
    ap.add_argument('--max-size', type=int, default=2000)
    ap.add_argument('--fdr', type=float, default=0.05)
    args=ap.parse_args()
    inp, out=Path(args.input).resolve(), Path(args.output).resolve()
    annot=find_annotation(inp)
    groups=discover_deg_pairs(inp)
    if not groups:
        raise FileNotFoundError(f'No DEG input files found in {inp}. Supply e.g. 1dpi.up.txt and 1dpi.down.txt.')
    bg=Path(args.background).resolve() if args.background else inp/'background.txt'
    if not bg.is_file():
        bg=None
        print('WARNING: background.txt is absent. Using eggNOG annotation genes as background. For published RNA-seq ORA use all statistically tested genes instead.',file=sys.stderr)
    # Read DEG tables (gene_id column, other recognized ID columns, or 1 ID/line).
    gene_sets={key:read_gene_table(path,None,False) for key,path in groups.items()}
    annot_ids = set(read_emapper(annot)['query'].astype(str).str.strip())
    input_ids = set().union(*gene_sets.values())
    id_regex, id_report = choose_id_rule(annot_ids, input_ids, args.id_mode, args.id_regex)
    print('Gene ID matching: ' + json.dumps(id_report), flush=True)
    if id_report.get('matched_after_rule', id_report.get('exact_matches', 0)) == 0:
        raise ValueError('No overlapping gene IDs between annotation and DEG files. Check identifiers or supply --id-mode regex --id-regex PATTERN.')
    gene_sets = {key: {normalize_gene(g, id_regex) for g in genes} for key, genes in gene_sets.items()}
    out.mkdir(parents=True,exist_ok=True)
    (out/'id_matching.json').write_text(json.dumps(id_report,indent=2),encoding='utf-8')
    # Construct all = unique(up U down) for EACH timepoint/comparison.
    # Write normalized lists for every direction so the same IDs are used for
    # database matching and enrichment, and the originals remain untouched.
    all_groups = {}
    derived=[]
    for comparison in sorted({name for name, _ in groups}):
        up=gene_sets.get((comparison, 'up'), set())
        down=gene_sets.get((comparison, 'down'), set())
        overlap=up & down
        if overlap:
            raise ValueError(f'{comparison}: {len(overlap)} genes appear in both up and down DEG lists, '
                             f'e.g. {sorted(overlap)[:5]}. Correct the input files before enrichment.')
        if not up or not down:
            print(f'NOTE: {comparison}: missing or empty direction; all will include the available genes.',file=sys.stderr)
        destination = out/'generated_DEGs'/comparison
        for direction, selected in [('up',up),('down',down),('all',up|down)]:
            file_path=destination/f'DEGs.{direction}.txt'
            write_gene_list(file_path, selected)
            # Empty groups are still written, but no enrichment can be run on them.
            if selected:
                all_groups[(comparison, direction)]=file_path
            derived.append(dict(comparison=comparison, direction=direction, unique_genes=len(selected),
                                generated_file=str(file_path),
                                source_up=str(groups.get((comparison,'up'),'')),
                                source_down=str(groups.get((comparison,'down'),''))))
    pd.DataFrame(derived).to_csv(out/'generated_deg_summary.tsv',sep='\t',index=False)
    # Resolve annotations automatically; cache official references in resources/.
    # Always use complete resources. Never silently rely on the example subset.
    go_path = Path(args.go_obo).resolve() if args.go_obo else ROOT/'resources'/'go-basic.obo'
    kegg_path = Path(args.kegg_names).resolve() if args.kegg_names else ROOT/'resources'/'kegg_pathway_names.tsv'
    if args.go_obo:
        if not go_path.is_file(): raise FileNotFoundError(go_path)
    else:
        get_resource(go_path,RESOURCE_URLS['go'],valid_go,'GO ontology',args.offline)
    if args.kegg_names:
        if not kegg_path.is_file(): raise FileNotFoundError(kegg_path)
    else:
        get_resource(kegg_path,RESOURCE_URLS['kegg'],valid_kegg,'KEGG pathway names',args.offline)
    resource_manifest(out,go_path,kegg_path)
    dbdir=out/'database'
    print(f'\n[1/4] Build database from {annot.name}',flush=True)
    build(SimpleNamespace(emapper=str(annot),outdir=str(dbdir),go_obo=str(go_path),kegg_names=str(kegg_path),
          id_regex=id_regex,include_go_roots=False,force=True))
    report=[]
    for i,((comparison,direction),path) in enumerate(sorted(all_groups.items()),start=2):
        print(f'\n[{i}/{len(all_groups)+1}] Enrichment: {comparison} / {direction}',flush=True)
        folder = out/direction if comparison == 'default' else out/comparison/direction
        enrich(SimpleNamespace(database=str(dbdir/'annotation.sqlite'),deg=str(path),background=str(bg) if bg else None,
              gene_column=None,background_column=None,padj_column=None,padj_cutoff=0.05,lfc_column=None,
              direction='all',universe=args.universe,min_size=args.min_size,max_size=args.max_size,
              min_overlap=1,fdr=args.fdr,outdir=str(folder)))
        summary=json.loads((folder/'enrichment_summary.json').read_text())
        for cat, metrics in summary['categories'].items():
            report.append({'comparison':comparison,'deg_group':direction,'database':cat,
                           'submitted_DEGs':summary['input_unique_degs'],
                           'matched_DEGs':metrics['deg_size'],'background_genes':metrics['background_size'],
                           'terms_tested':metrics['terms_tested'],'significant_terms':metrics['significant_terms']})
    pd.DataFrame(report).to_csv(out/'enrichment_overview.tsv',sep='\t',index=False)
    print(f'\nFinished. Browse {out}/. Overview: {out/"enrichment_overview.tsv"}',flush=True)

if __name__=='__main__':
    main()
