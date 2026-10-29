# DuoGran TPDS manuscript

LaTeX source for the DuoGran paper in IEEE Computer Society journal format.

## Build

Install a LaTeX distribution with IEEEtran (including its bibliography style), pdfLaTeX, BibTeX, and latexmk, then run:

```bash
make
```

The generated manuscript is `Source/prolet.pdf`. It is not tracked in Git.

The repository includes the main manuscript (`Source/prolet.tex`), bibliography (`Source/reference.bib`), build files, and all figure PDFs referenced by the manuscript. The figure PDFs are build inputs and are tracked.

## Figure 5 generation

Plotting code and derived data are managed in the separate prolet code repository:

- `prolet/analysis/plot_workload_fct.py`
- `prolet/analysis/figure_data/fig5_workload_fct.json`
- `prolet/analysis/README_paper_figures.md`

This manuscript repository keeps the generated `Source/figures/workload_fct_comparison.pdf` and the article text. It does not keep a duplicate plot generator or statistics snapshot.

With NumPy and Matplotlib installed, run from the prolet checkout:

```bash
python3 analysis/plot_workload_fct.py \
  --archive ../prolet-experiment-logs \
  --output output/paper-fig5
cp output/paper-fig5/workload_fct_comparison.pdf /path/to/TPDS/Source/figures/
```

The generator verifies raw-input/FCT hashes, paired baseline flow identities and ideal FCTs, and complete per-scheme size-bin coverage. Its output directory contains the PDF, PNG preview and regenerated `plot_data.json`. Updating the manuscript figure is an explicit copy; the script does not modify this repository automatically.

ConWeave and CONGA use complete paired 100-ms inputs. DuoGran retains the original 10-ms curves; flow cohorts, hashes, axis transformations and clipped points are recorded in the prolet statistics snapshot. Both workloads use a linear lower range and progressive upper-range compression, with all ConWeave curves fully visible.
