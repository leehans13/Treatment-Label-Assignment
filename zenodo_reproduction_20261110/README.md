# Reproducing the treatment-label assignment analyses

This package contains the analysis code for the manuscript's main figures. It is based on the DEGpaper analysis project: the Kang and Schmidt scripts use its shared model implementations, and the null/power scripts adapt its simulation notebooks to read the public H5AD files. The package produces model results and numerical summaries; it does not contain figure-drawing code or manuscript images.

The simulation scripts use fixed seeds to generate samples and record their results.

## Which analysis goes with each figure?

| Main figure | Analysis and output |
| --- | --- |
| Figure 1 | A conceptual summary of the study designs and method performance. It synthesizes the theory and the simulation and perturbation analyses below; there is no separate model fit that produces the table. |
| Figure 2 | Reed subject-level assignment simulation: `src/nullpower/prepare.R`, `run.R`, and `summarize.R` with dataset `bca`. The full analysis writes `outputs/null_subject/` and `outputs/pow_subject/`; the main figure uses one cell type. |
| Figure 3 | Reed cell-level assignment simulation: the same scripts and dataset, using `outputs/null_cell/` and `outputs/pow_cell/`. |
| Figure 4 | Kang interferon stimulation: `src/kang/run.R`; Schmidt CRISPR perturbation: `src/schmidt/run.R`. Each writes model results and pairwise rejection and Spearman matrices. |

The package also includes the Yazar and Reichart versions of the null/power simulations, the synthetic variance-scale analysis (`src/variance_scale/`), and the Schmidt non-targeting-control analysis (`src/schmidt/run_null_control.R`). These are supplementary analyses.

## Environment

Run commands from this package directory. Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), then create the environment and install NEBULA:

```sh
micromamba create -f environment.yml
micromamba run -n r-bio Rscript scripts/install-extra-packages.R
```

The source code has been smoke-tested on macOS arm64 with an existing `r-bio` environment. NEBULA requires a working C++ compiler. If an environment named `r-bio` already exists, use a different environment name in `environment.yml` rather than replacing it, and substitute that name in the commands below.

## Get the data

The input datasets are not included. The download script retrieves the Reed, Yazar, and Reichart H5AD files, the Schmidt archive, and two Kang Figshare count files. It resumes partial downloads and checks file sizes; it also checks the published MD5 hashes for the Schmidt and Kang files. Total download size is about 23 GB; extraction and analyses need additional space.

```sh
bash scripts/download-data.sh
```

The pinned sources are the [Reed](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637), [Yazar](https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1), and [Reichart](https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb) CELLxGENE collections; the [Schmidt CRISPRa dataset](https://zenodo.org/records/5784651); and the [Kang Figshare matrices](https://figshare.com/articles/dataset/PBMC_scRNAseq_Kang_2018_RDS_dgCMatrix/22572694). The exact file URLs, expected sizes, and available hashes are in the download script. The Figshare matrices do not include the annotations needed by the Kang analysis. Fetch the annotated ExperimentHub object separately:

```sh
mkdir -p data/kang
micromamba run -n r-bio Rscript -e 'saveRDS(muscData::Kang18_8vs8(), "data/kang/Kang18_8vs8.rds")'
```

For Schmidt, extract the two files used by the analysis:

```sh
unzip -n data/schmidt/CRISPRa-Perturb-seq.zip \
  'Perturb-seq/data/HuTcellsCRISPRaPerturbSeq_Re-stimulated.h5Seurat' \
  'Perturb-seq/data/GO_0005125_Cytokines.txt' -d data/schmidt
```

The simulation reader uses Reed and Reichart `raw/X` counts and Yazar `X`. The Kang analysis uses the annotated ExperimentHub object. The Schmidt analysis selects normalization genes using SCT counts across all cells.

## Run the analyses

The full jobs are computationally expensive. They write under `data/` and use checkpoints; rerun the same command to resume, and use a fresh output directory after changing inputs, code, or software versions. Paths and modification times are part of checkpoint validation.

### Null and power simulations

Prepare each downloaded dataset once, then run the full analysis and summarize it. Replace `bca` with `yazar` or `reichart` for the other datasets:

```sh
micromamba run -n r-bio Rscript src/nullpower/prepare.R bca
micromamba run -n r-bio Rscript src/nullpower/run.R bca full
micromamba run -n r-bio Rscript src/nullpower/summarize.R bca
```

The full analysis runs six cell types, sample sizes 10 and 50, both label-assignment scenarios, and ten iterations: 240 checkpoints per dataset. Each iteration uses 100 genes, including five altered genes. The summary writes per-condition p-value CSVs, rejection summaries and model-error reports under `data/nullpower-reproduction-20261104/bca/outputs/`. It refuses to overwrite a nonempty output directory.

A small functional check, after preparing the input, uses one cell type, ten donor draws, 20 genes and one iteration:

```sh
micromamba run -n r-bio Rscript src/nullpower/run.R bca smoke
```

A cell-assignment check is available as `smoke-cell`. Smoke results are stored separately from the full run.

Donors are sampled **with replacement**; repeated draws receive separate analysis IDs and can reuse cells. The same selected genes are used across iterations for a setting.

### Perturbation data

```sh
micromamba run -n r-bio Rscript src/kang/run.R
micromamba run -n r-bio Rscript src/schmidt/run.R
```

Kang writes `kang.result.rds`, `kang.pow.csv`, and `kang.sp.csv` under `data/kang/results/`. Schmidt writes `schdmidt.result.rds`, `schmidt.pow.csv`, and `schmidt.sp.csv` under `data/schmidt/results/`. The CSV matrices summarize shared rejections and Spearman correlations using the notebook's thresholds and non-finite-value handling. The scripts accept optional input and output paths as shown in their first lines. The notebook analyses use NEBULA, three pseudobulk methods, and two cell-wise methods; the Schmidt run uses one NEBULA worker. To check each analysis on a small real-data subset without running the full job, use `micromamba run -n r-bio Rscript scripts/smoke-perturb.R kang` or `micromamba run -n r-bio Rscript scripts/smoke-perturb.R schmidt`.

### Other analyses

```sh
micromamba run -n r-bio Rscript src/variance_scale/run.R
micromamba run -n r-bio Rscript src/schmidt/run_null_control.R
micromamba run -n r-bio Rscript src/schmidt/summarize_null_control.R
```

The variance simulation defaults to 100 replicates per setting and 2,000 bootstrap draws. A quick executable check is `Rscript src/variance_scale/run.R /tmp/variance-smoke 2 10 smoke 1` in the environment. The Schmidt control analysis defaults to 69 random splits of the same control cells. Optional arguments are described at the top of the scripts.
