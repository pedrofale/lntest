# Reproducing the results

Everything in the Nature submission that depends on code lives here, one directory per experiment.
This file is the reproduction instruction; the top-level `README.md` covers the `lntest` package itself.

Read [Known gaps](#known-gaps) before trusting any number you produce.

## Run everything from this directory

```bash
cd reproducibility
python -m <arm>.<script>          # e.g. python -m synthetic_nb.small_test
```

Not `python synthetic_nb/small_test.py`.
The `-m` form puts the **working directory** on `sys.path`, so a script inside an arm can import the shared modules that sit beside the arms — `paths`, `baselines`, `de_utils`, `evaluation_utils`, `utils_frozen` — with no `__init__.py`, no packaging file for `reproducibility/`, and no `sys.path` manipulation.
Running a script by path puts *its own directory* on `sys.path` instead, and the shared imports fail.

**Notebooks are the exception.** Run a notebook from the arm directory that contains it; each does `sys.path.append('../')` in its first cell to reach the same shared modules.

## The environment

```bash
conda env create -f env.yaml
conda activate de-ziln-reproducibility
pip install -e ..                 # installs the lntest package from ../src/
export R_HOME="$CONDA_PREFIX/lib/R"
```

**`R_HOME` is not optional.** Without it, `import rpy2.robjects` aborts the interpreter with

```
Error in substring(x, m + 1L) : invalid substring arguments
```

and exit code 139, *before any of your code runs*, because `clustering/de.py`, `de_utils.py` and `lymphnode/subsampling.py` all import rpy2 at module scope.
`--skip-mast` does not help: the import happens either way.
This affects the clustering, celltype and lymph node arms. R itself is fine — only rpy2's initialisation of it needs the variable.

**Figures are written as PDF.** Every plotting entry point in this tree routes its filename
through `paths.fig_name`, which forces the extension. Set `LNTEST_FIG_FORMAT` to override for a
run (`LNTEST_FIG_FORMAT=png python -m theory.concave_ordering`).

The tree used to emit a mix of `.eps` and `.png`, and four scripts hardcoded `format='eps'` in the
`savefig` call, so `--output x.png` produced PostScript named `.png`. Those are gone; the filename
now decides. PDF is also markedly smaller for these plots — `variance_vs_fpr_mu5` went from 460 KB
raster to 17 KB vector.

**Two consequences worth knowing before regenerating anything.** The 22 figure entries in
`reference/clustering/checksums.sha256` name `.png` files and cannot match PDF output; they were
already due for regeneration after the single-trigamma change, and 10 of them were never
verifiable anyway. And any `\includegraphics` in the manuscripts that spells out `.png` or `.eps`
needs updating, or the extension dropped so LaTeX picks what exists.

**`pandas<3` is a hard pin, not caution.** pandas 3.0 makes the arrow-backed `str` dtype the
default for string columns and indices, and anndata 0.12.6 has no writer registered for it, so
`adata.write()` fails on *any* AnnData:

```python
import numpy as np, anndata as ad
ad.AnnData(np.ones((3, 2))).write("probe.h5ad")   # IORegistryError
```

Reading is unaffected, so an arm loads its input, computes everything, and dies on the last line.
`celltype/de.py` does exactly that: it finishes LN, t-test, Wilcoxon and MAST, then fails on its save.
`clustering/de.py` is the one arm that survives pandas 3, because it carries its own
`coerce_arrow_strings()` helper (added in `c7c5d94`, the commit behind the reference checksums) and
calls it immediately before writing. Do not read that as the problem being solved — it is solved in
one file. If you built this environment before the pin was added,
`conda install -n de-ziln-reproducibility 'pandas<3'`.

**If your environment predates 2026-09-10**, update it: `geopandas`, `shapely` and `pyarrow` were missing, so neither spatial arm's preprocessing could run.

```bash
conda env update -f env.yaml --prune
```

## What you can run right now

Three data files are committed for the live arms — `citeseq/data/memory_CD4.h5ad` and the two
`nullsplit/data/split_*.csv` — plus five under `synthetic_nb/legacy/data/`. Everything else is
fetched, downloaded or generated.

| Arm | Input | In a fresh clone? |
|---|---|---|
| `synthetic_nb` | none — generates its own | **yes** |
| `theory` | none | **yes** |
| `nullsplit` | `nullsplit/data/split_{a,b}.csv` (committed) | **yes** |
| `citeseq` | `citeseq/data/memory_CD4.h5ad` (committed, 16 MB) | **yes** — also rebuildable from a download, see below |
| `clustering` | `data/pbmc3k_filtered_gene_bc_matrices.tar.gz` | **yes** — `python -m fetch_data --arm clustering` |
| `celltype` | `data/pbmcs3k_pre.h5ad`, `data/kang_2018.h5ad` | **yes** — `fetch_data --arm celltype`, then `python -m celltype.preprocess` |
| `kidney` | `kidney/data/{merged_blobs_in_cluster_5,podocytes_2um}.h5ad` | **download yes** — `fetch_data --arm kidney` (6.5 GB), then `preprocess.ipynb` |
| `lymphnode` | `lymphnode/data/vishd-cluster1-cluster3-2um-with-clusters-subsampled.h5ad` | **download yes** — `fetch_data --arm lymphnode` (4.7 GB), then the reproduce script (hours) |

`data/` and `*/data/` are gitignored, so a file being present in one working copy says nothing about a fresh clone.
Every entry point checks its input first and exits naming the missing path and where it comes from, rather than failing part-way through.

Fetching them:

```bash
python -m fetch_data                  # everything missing, verified against sha256
python -m fetch_data --arm clustering # just what one arm needs
python -m fetch_data --list           # what is present, what is not
```

`data_sources.yaml` carries a URL for every fetched input — all seven, spatial included — and
`fetch_data.py` downloads what is missing, verifies it, and unpacks the tarballs.

Verification is two-tier. Four entries are pinned by **sha256**. The three multi-GB spatial archives
carry a **byte count only**, because computing a digest for them means downloading 11 GB first; those
report as `SIZE-OK` rather than `OK`, which catches truncation and nothing else. Pin them properly by
hashing after a download and filling in the `sha256:` field. It is idempotent, and it
refuses to overwrite a file whose digest does not match rather than silently replacing it.

**CITE-seq rebuilds from its download too**, as of 2026-09-18:

```bash
python -m fetch_data --arm citeseq          # 30 MB raw 5' PBMC 10k matrix
cd citeseq/R && Rscript pbmc10k_process.R   # ~17 s -> results/pbmc10k_cd4_memory.rds
Rscript pbmc10k_to_h5ad.R                   # ~2 min -> data/memory_CD4.h5ad
```

Verified: the rebuilt `memory_CD4.h5ad` matches the committed one exactly — 1318 x 2035, same cells,
same genes, `max|dX| = 0`. That also settled which raw file is the source: 10x serves two under the
same name, and only `cell-vdj/5.0.0` reproduces it, so the entry is now sha256-pinned.

`memory_CD4.h5ad` stays committed for convenience — the rebuild takes ~2 minutes and needs the R
stack — but it is no longer the only way to get it. Note the rebuilt file is *content*-identical and
not *byte*-identical: HDF5 containers differ in compression and metadata, so `git status` will show
it as modified after a rebuild even when nothing in the data changed.

One thing this does **not** do: **no Zenodo DOI exists**, so the non-derivable intermediates named in
the vault's `decision-data-out-of-head` have nowhere to be fetched from.

Verified 2026-09-18: deleting `data/pbmc3k*` and re-running `python -m fetch_data --arm clustering`
downloads, verifies and unpacks it, after which `clustering.de` and `celltype.preprocess` both run
off the fetched matrix.

`data/pbmcs3k_pre.h5ad` is the pbmc3k matrix after standard QC; `clustering/legacy/resolution_sweep.ipynb` writes it, and `clustering/de.py` applies the same QC inline when handed the raw `.tar.gz`.

## The arms

### `synthetic_nb` — negative-binomial simulations

```bash
python -m synthetic_nb.small_test      # ~10 s. The reviewer-facing single run
python -m synthetic_nb.de_test         # -> synthetic_nb/results/nde_mu10_be/
python -m synthetic_nb.latex_tables    # formats de_test's CSV
python -m synthetic_nb.null
python -m synthetic_nb.lfc_confidence_intervals
```

`small_test.py` reproduces **one** run of a table averaged over 20; its numbers are not expected to match the paper exactly.
The NB parameters behind the paper's table are not captured in any config: the submission said they "need to be adjusted according to the text", and nothing here records what they were.

### `citeseq` — CITE-seq, surface protein as ground truth

```bash
python -m citeseq.exp                  # writes to output/citeseq/
```

The R stages that build `memory_CD4.h5ad` from the raw 10x download are in `citeseq/R/`, run with `Rscript`.
The committed `memory_CD4.h5ad` means you do not need them unless you are rebuilding from raw.
`citeseq/R/pbmc10k_seurat.R` generates Seurat LFC estimates that **were not used in the paper**, because they did not affect the conclusion; it is kept for provenance.

### `clustering` — PBMC3k, clustering resolution sweep

```bash
python -m clustering.de      --config clustering/config.yaml --output-dir output/clustering

# metrics and plots do NOT read --output-dir to find their input. Point them at
# what the previous stage wrote, or they fail / silently do nothing. See below.
sed 's|^adata_path:.*|adata_path: "output/clustering/pbmc3k_filtered_gene_bc_matrices_DE.h5ad"|' \
    clustering/config.yaml > .metrics_config.yaml          # gitignored, see .gitignore
python -m clustering.metrics --config .metrics_config.yaml --output-dir output/clustering
python -m clustering.plots   --metrics-dir output/clustering/metrics --config clustering/config.yaml
```

Run them in that order, and **pass `--skip-mast`** (it still needs `R_HOME` — the flag does not avoid
the rpy2 import). Measured 2026-09-10:

| | DE stage | Reference metric CSVs | Reference figures |
|---|---|---|---|
| with MAST | **2 h** (63 fits, one per cluster; the 30-cluster resolution is most of it) | 68/68 byte-identical | 12/22 |
| `--skip-mast` | **15 s** | 68/68 byte-identical | 12/22 |

The two are **equivalent for every verifiable output**, because nothing downstream reads the MAST
results. `metrics.py` hardcodes `de_methods = ["ln", "t_test", "wilcoxon"]` and neither it nor
`plots.py` mentions MAST anywhere, so the 5 `mast_*` keys the 2-hour stage writes into the `_DE.h5ad`
are stored and never consumed. Run with MAST only if you need those keys for something outside this
tree — they are not deleted here, because whether the manuscript leans on them is an open question.

This is the arm with reference checksums: see [Outputs](#outputs-and-reference-values).
Verified end to end on 2026-09-10 — **all 68 metric CSVs byte-identical to `reference/`**.

Two wiring defects to know about, both found by running it:

- **`metrics` re-reads `config['adata_path']` and opens it with `read_h5ad`.** It does not read the
  `_DE.h5ad` that `de` just wrote into `--output-dir`. With the committed config, which names the
  `.tar.gz`, it dies with `OSError: file signature not found`. Hence the `sed` above.
- **`clustering/plots.py` silently drops four figures without `--config`.** Line 75 reads
  `top_genes = config.get('top_genes', 10)` while `clustering/config.yaml` sets `top_genes: 20`,
  so with no config the script looks for `avg_jaccard_top10_by_resolution.csv`, finds nothing, and
  its `boxplot()` helper returns early. You get 18 figures instead of 22, with no message. The
  heatmaps are unaffected because they `glob` rather than build the name, which is why the `top20`
  heatmaps appear while the `top20` boxplots do not. One of the four,
  `avg_jaccard_top20_by_resolution`, is a panel in the Nature submission. Always pass
  `--config clustering/config.yaml`.
- **`plots --metrics-dir` wants the `metrics` directory itself**, not the output root — its `--help`
  says `e.g. output/metrics`. Given the output root it finds no CSVs, prints
  `(No recall_matrix_*.csv files found; skipping ...)`, reports `Plots complete!`, **exits 0 and
  writes nothing**. Figures default to a `figures/` sibling of `--metrics-dir`, which is where the
  reference manifest expects them.

**Ten of the 22 figures cannot be byte-reproduced by anyone**, including whoever made the reference.
The `*_by_resolution.png` boxplots overlay `sns.stripplot`, whose jitter is drawn from the unseeded
global RNG, so two consecutive runs of the same script on identical CSVs differ. Confirmed by running
it twice: the 12 heatmaps are byte-identical both to each other and to `reference/`; the same 10
boxplots differ every time. Their reference hashes are therefore not verifiable — judge those panels
by the metric CSVs underneath them, which do match exactly.

### `celltype` — Kang 2018 cell types

```bash
python -m celltype.de --config celltype/config.yaml --output-dir output/celltype
```

`celltype/kang_markers.ipynb` produces the committed marker CSVs in `celltype/results/`.

### `kidney` — Visium HD, glomerular capsules

Preprocessing first, either form -- they do the same work:

```bash
python -m kidney.preprocess      # from reproducibility/
```

or open `kidney/preprocess.ipynb` **from `reproducibility/kidney/`** and run it top to bottom.
Either downloads 6.5 GB from 10x and writes both matrices into `kidney/data/`.
Took ~9 min on an M-series laptop; produces 161 capsules over 18,085 genes, which is what
`main_recomb25.tex` reports for this dataset.

```bash
python -m kidney.umi_null   --n_cells_remove 35 --n_reps 10000
python -m kidney.umi_de     --n_cells_remove 35 --q 0.005 --lfc 3 --p_min 0.125 --n_reps 10000
python -m kidney.spot_split --n_shape_ids_remove <unrecorded>
python -m kidney.fpr_plots        --csv_file <umi_null results.csv>
python -m kidney.umi_de_plots     --input <umi_de results.csv> --output de.eps
python -m kidney.spot_split_plots --csv_file <csv> --json_file <json>
```

**These parameters come from the manuscript, not from the code.** `--n_cells_remove`, `--q` and
`--lfc` are `required=True` with no defaults, so nothing in this repository records what was run;
until 2026-09-24 the block above showed `50 / 0.1 / 1.0`, which were placeholders and matched no
published figure. The values now shown are what `main_recomb25.tex` states:

| Flag | Value | Source |
|---|---|---|
| `--n_cells_remove` | 35 | "we first filtered out the 35 capsules with the lowest total UMI" (§ Visium HD); the panels are also named `fdr_results_35.eps` and `DE_plot_lfc3.eps` |
| `--q` | 0.005 | "With probability 0.005, we sample DEGs ... around 90 DEGs" |
| `--lfc` | 3 | "amplified to ensure a three-fold $\log_2$ change, i.e. LFC=3" |
| `--n_reps` | 10000 | "For each downsampling ratio $p$, 10,000 random tests were ran" |
| `--p_min` | 0.125 | "$p\ge 2^{-3}$" -- a **derived** bound, see below |

`p_min` is not free: the manuscript ties it to `lfc` as $p \ge 2^{-\mathrm{LFC}}$, because a gene
amplified by LFC log2-units cannot survive downsampling below that ratio. **No code enforces this.**
Change `--lfc` and you must change `--p_min` with it; the script defaults (`p_min 0.1`) satisfy the
constraint for LFC=3 only by being close to 0.125, and violate it for any LFC below 3.

`spot_split`'s capsule filter is left unfilled on purpose: the spot-subsampling section states no
such filtering, so there is no published value to copy. Its `p_min 0.1 / p_max 0.5` defaults do
match the paper's illustrated splitting probabilities of 0.5, 0.3 and 0.1.

**Full-depth check, 2026-09-24.** Before any subsampling, splitting all 161 capsules by whether
total UMI exceeds the median library size reproduces the manuscript's negative-control counts:

| Method | Here | `main_recomb25.tex` |
|---|---|---|
| LN's $t$-test | 10 | 10 |
| log1p $t$-test | 1099 | 1099 |
| Wilcoxon | 2878 | 2876 |

Two exact, Wilcoxon +2 of 2876 (0.07%), most likely rank-tie handling across scipy versions.
Takes seconds, and is worth running before committing to a 10,000-replicate sweep.

The two UMI scripts read `merged_blobs_in_cluster_5.h5ad` (one row per capsule); `spot_split` reads `podocytes_2um.h5ad` (one row per 2 µm spot, with a `shape_id` column). They are not interchangeable — see [Known gaps](#known-gaps).

### `lymphnode` — Visium HD, Cluster-1 vs Cluster-3

Preprocessing is a driver script; read its header first, it downloads 4.4 GB and the middle step is hours of single-threaded point-in-polygon:

```bash
bash lymphnode/reproduce_vishd_cluster1_cluster3.sh
```

Then:

```bash
python -m lymphnode.subsampling --config lymphnode/config_50rep.yaml --output-dir output/lymphnode
python -m lymphnode.ln_de_vs_rest   --input <h5ad> --output de.csv
python -m lymphnode.cluster_de_gsea --input <h5ad> --output_de de.csv --output_gsea gsea.csv
```

### `theory` — illustrations, no data

```bash
python -m theory.mean_ci_coverage    # ~1 s
python -m theory.lfc_ci_coverage
python -m theory.concave_ordering    # ~1 s
```

These are exempt from depending on `lntest` — being readable in one file, with the algebra inline, matters more here.
The exemption is on the *code*, not the maths: where the manuscript prints a formula, that formula is the specification.

### `nullsplit` — null split false-positive rate

`nullsplit/ziln_null_fpr.ipynb`, run from `reproducibility/nullsplit/`. Its two CSVs are committed.

## Outputs and reference values

Results belong in `output/`, which is gitignored: `output/clustering/`, `output/lymphnode/` and
`output/citeseq/`. The synthetic NB and theory arms write to their own `results/`, whose outputs git ignores.
**Figures, tables and checksums are not committed** — `.gitignore` excludes `*.pdf`, `*.csv`, `*.tex`
and `*.sha256`; files that were already tracked stay tracked.

**Every output is byte-reproducible except wall-clock timings.** Regenerated twice from scratch on
2026-09-24, under Bonferroni for every gene-wise test and LN's trigamma as psi_1(a) = 1/a, and checked
file by file. Three things make that hold, and each is easy to undo by accident:

- PDFs carry a timestamp unless pinned. `paths.py` sets `SOURCE_DATE_EPOCH` for matplotlib; the R
  scripts rewrite the date after `ggsave`, because R's PDF device ignores that variable.
- `clustering/plots.py` seeds the global RNG before each stripplot, which is where seaborn draws its jitter.
- `lymphnode/subsampling.py` sorts its results before writing, because workers finish in any order.

The exceptions: `lymphnode/de_time_s_vs_subsampling.pdf`, the `de_time_s` columns of the lymph node CSVs,
and a timing field in its metadata JSON.

The checksums of the current outputs live in `output/reference/`, outside git; verify one with
`shasum -a 256 -c output/reference/<arm>/checksums.sha256`. The committed `reference/` holds the
June 2026 manifests, which pin the run behind the published PBMC3k and CITE-seq figures.

| Manifest | Files |
|---|---|
| `output/reference/clustering/` | 68 metric CSVs, 22 figures |
| `output/reference/citeseq/` | metrics, per-gene results, LaTeX table, 102 figures |
| `output/reference/lymphnode/` | 28 figures (the CSVs carry wall-clock columns, so compare them by value) |
| `output/reference/synthetic_nb/` | the two variance sweeps and three dispersion boxplots |
| `output/reference/theory/` | two figures, two CSVs |

`reference/unattributed/` holds two `.npy` files no script reads. Read `PROVENANCE.md` before assuming anything about them.

The regenerated manuscript panels are staged under `output/manuscript_figures/`, with Overleaf paths
and a README of what each one is. Reaching the paper is a manual upload.
The BH-era outputs they replace are kept in `output/archive/`.

**Re-running an arm can overwrite committed files.** Nine tracked files live under `*/results/` rather than `output/`:

```
celltype/results/kang_{B,T}_markers.csv
citeseq/results/CITE_seq_{results,metrics}.csv
citeseq/results/CITE_seq_lfc_error_plot.png
citeseq/results/CD{3,4,45RA}_density_plot.pdf
synthetic_nb/results/de_metrics_mu10_nobatch.csv
```

The CITE-seq arm used to write straight over its tracked CSVs; it now writes to `output/citeseq/`, and the tracked files keep the June outputs. Check `git status` after any run all the same.

## The equivalence harness

`check_utils_vs_lntest.py` compares the frozen RECOMB-era estimator in `utils_frozen.py` against the published `lntest` package, arm by arm.

```bash
python -m check_utils_vs_lntest --arm all
```

It exists because the arms were migrated onto `lntest` one at a time and something had to prove no headline number moved.
`utils_frozen.py` is a frozen copy, kept for that comparison and for the two `theory/` scripts — do not build on it.

## Known gaps

- **Nothing here has a test suite.** `small_test.py` is a reviewer-facing single run, explicitly not expected to match the submission.
- **`theory/concave_ordering.py`, `theory/lfc_ci_coverage.py` and `synthetic_nb/null.py` call `plt.show()` and save nothing.** Run them in a notebook, or add a `savefig`; headless they complete and produce no file.
- **The kidney arm's subsampling parameters were never recorded in code.** `--n_cells_remove`, `--q` and `--lfc` are required arguments with no defaults and no config, so the only record of what produced the published panels is the manuscript prose. Recovered 2026-09-24 and written into the command block above; `spot_split`'s capsule filter has no stated value anywhere.
- **The kidney arm's published spot-split panel has unclear provenance.** `spot_split.py` defaulted to the per-capsule aggregate, whose `shape_id` is an index rather than a column, so its own guard raised `ValueError` before doing any work. The default now points at `podocytes_2um.h5ad`, per this repository's own top-level README, but which file produced the published figure is an open question.
- **`celltype/config.yaml` names an input that cannot satisfy it.** It asks for
  `celltype_column: cell_type` with the PBMC3k labels `Naive CD4 T` and `CD14 Mono`, in
  `data/pbmcs3k_pre.h5ad` — whose `obs` is only `n_genes, percent_mito, n_counts`. Nothing in this
  tree writes a `cell_type` column onto that file, and the config has not changed since it was first
  committed, so the arm has never run as configured. `de.py` itself is fine: given `data/kang_2018.h5ad`
  and two of its real labels it runs LN, t-test, Wilcoxon and MAST and saves. Which input was meant is
  an open question.
- **The null split notebook does not reproduce its committed outputs.** `ziln_null_fpr.ipynb` shows 10
  rejections / FPR 0.00304; running it gives 21 / 0.00638. The notebook was last executed 2024-11-06 and
  the estimator changed nineteen times after that, so its outputs are stale rather than wrong — 21 is what
  the frozen RECOMB-era estimator gives. See question 9 of the vault's math-questions handoff.
- **Neither spatial arm has been run end to end here.** The paths are consistent and the inputs are public, but nobody has spent the 11 GB.
- **The trigamma question is settled.** `lntest` has one trigamma difference, `trigamma_diff(a, n) = 1/a - 1/n`, and no flag to select another. That is the form the paper defines — `psi_1(z) = 1/z`, at `nature_submission/sections/methods.tex:79-83` — and the form every published result came from. Nothing in this tree asks for a trigamma any more, because there is nothing to ask for.
