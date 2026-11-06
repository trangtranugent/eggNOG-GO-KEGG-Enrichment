# eggNOG GO & KEGG enrichment

**Author:** Trang Tran

A portable Python workflow for GO and KEGG pathway over-representation analysis (ORA) using an existing eggNOG-mapper annotation file and gene lists from RNA-seq. It works for any organism with usable eggNOG annotations; no organism-specific reference package is required.

## Quick start

```bash
python -m pip install -r requirements.txt
python run.py
```

The first run downloads the full GO ontology and KEGG pathway-name list to `resources/`; future runs reuse these resources. **No separate resource commands are needed.** An Internet connection is required for the first run. For reproducible reruns, preserve the cached files and `output/resource_manifest.json`.

Put an annotation file and **any** up/down DEG text files in `input/`:

```text
input/
├── emapper.annotations.tsv     # or emapper.annotations.xlsx (only one)
├── background.txt              # all RNA-seq-tested genes (strongly recommended)
├── DEGs.up.txt                 # simple one-comparison case
├── DEGs.down.txt
├── 1dpi.up.txt                 # arbitrary comparison/timepoint names
├── 1dpi.down.txt
├── TreatmentA_3dpi.up.txt
└── TreatmentA_3dpi.down.txt
```

**No `DEGs.all.txt` is needed.** The program finds every `*.up.txt` and `*.down.txt`, pairs files with the same prefix, and generates `all = unique(up ∪ down)` for each comparison. Prefixes `DEGs.1dpi` and `1dpi` both work; the `DEGs.` prefix is optional. Input files are not modified. Missing directions are allowed (but an empty direction is not enriched). Genes occurring in both directions for one comparison cause an error, to prevent contradictory results. Existing `*.all.txt` files are ignored.

Generated lists are saved under `output/generated_DEGs/<comparison>/`, with `DEGs.all.txt`, `DEGs.up.txt`, and `DEGs.down.txt`. Enrichment results are saved under `output/<comparison>/<direction>/`, or `output/<direction>/` for plain `DEGs.up.txt`/`DEGs.down.txt`.

The eggNOG file must have a `query` column and ideally `GOs`, `KEGG_Pathway`, and `KEGG_ko` columns. The `.txt` gene files can be one ID per line or TSV with a `gene_id` column. **DEGs.all.txt is not the statistical background.** Define DEG significance thresholds upstream (e.g., adjusted *p* value and log2FC). When `background.txt` is absent, annotation genes are used with a warning; provide a tested-gene background for publishable RNA-seq results.

## Run the supplied examples

```bash
python run.py --input example/auto_timecourse --output output/example_run
python -m unittest discover -s tests -v
```

The `example/auto_timecourse/` folder contains `1dpi.up.txt`, `1dpi.down.txt`, `3dpi.up.txt`, and `3dpi.down.txt` (no pre-made all files), with a synthetic eggNOG annotation and tested-gene background. It is ready to run.

Both examples are **synthetic** and show the expected input format, not biological findings.

## Gene ID matching (works across organisms)

The default `--id-mode auto` compares annotation and DEG identifiers. It **preserves full IDs** unless a recognizable terminal transcript suffix can be removed to increase matches. Common patterns include `_T001`, `.t1`, `-T1`, `.p1`, and `_i1`. The same transformation is applied to annotation, DEG and background identifiers. The chosen rule and matching counts are written to `output/id_matching.json`.

```bash
python run.py --id-mode exact   # always preserve all ID suffixes
python run.py --id-mode regex --id-regex '_isoform[0-9]+$' # custom pattern
```

There is **no truly universal transcript suffix**; if the annotation IDs and RNA-seq IDs use unrelated namespaces (e.g. protein accessions vs gene IDs), a separate gene ID mapping is needed. Check the match counts and avoid overly broad regexes. `auto` will not invent a mapping when no candidate improves matches.

## Results

```text
output/
├── id_matching.json
├── database/                  # SQLite + gene2go, gene2kegg, gene2ko mappings
├── generated_DEGs/             # automatically generated up/down/all gene lists
├── 1dpi/
│   ├── all/                   # go/kegg_all_results.tsv and *_significant.tsv
│   ├── up/
│   └── down/
├── 3dpi/
└── enrichment_overview.tsv
```

GO and KEGG analyses use one-sided hypergeometric enrichment with Benjamini–Hochberg FDR correction, separately for each annotation category. Output includes term IDs, overlap counts and gene IDs, adjusted *p* values, and gene background sizes. **This is ORA, not ranked GSEA.** GO names and namespaces, and KEGG pathway names, are retrieved automatically on first run from the official GO ontology and KEGG REST reference list. These are cached in `resources/` and reused on subsequent runs. A `resource_manifest.json` records each file’s source and SHA-256 checksum. Full reference files are not redistributed in this GitHub package; Internet is required for the first run (or supply local reference files using `--go-obo` and `--kegg-names`). Use `--offline` to forbid network access once the cache exists. This affects descriptions only, not enrichment statistics. See `python run.py --help` for options.

## References

- Cantalapiedra CP et al. (2021). *eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale*. Molecular Biology and Evolution 38:5825–5829. [doi:10.1093/molbev/msab293](https://doi.org/10.1093/molbev/msab293).
- Benjamini Y, Hochberg Y (1995). *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing*. Journal of the Royal Statistical Society B 57:289–300. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
- Please also cite the appropriate eggNOG **database version** and annotation tool version used to generate your input (which may differ from this workflow).

## Reference data and attribution

- GO ontology (`go-basic.obo`): [Gene Ontology Consortium](https://geneontology.github.io/docs/download-ontology/), retrieved from the official PURL. Record its version in publications.
- KEGG pathway labels: [KEGG REST API](https://rest.kegg.jp/), `list/pathway`. KEGG data are subject to the [KEGG terms/licensing](https://www.kegg.jp/kegg/legal.html); this repository downloads reference data for users rather than redistributing a bundled KEGG database.
- Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics*. [doi:10.1093/genetics/iyad031](https://doi.org/10.1093/genetics/iyad031).
- Kanehisa M et al. (2023). KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Research*. [doi:10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963).

## License

MIT © 2026 Trang Tran. See [LICENSE](LICENSE).
