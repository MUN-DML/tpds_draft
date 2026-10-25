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

The generator is `research/figures/fig5_workload_comparison/plot_fig5.py`; `plot_data.json` in the same directory records the plotted statistics, input hashes, bin definitions, axis transformation, and clipped points. It is a statistics snapshot, not an input consumed by the generator.

With NumPy and Matplotlib installed, run:

```bash
python3 research/figures/fig5_workload_comparison/plot_fig5.py
```

The script reads the local simulation results identified by `RUNS` and `WORKLOADS`; adjust `RUNS` for another checkout. The raw `run_out_fct.txt` files are not included. It verifies matching flow identities and ideal FCTs, regenerates the statistics snapshot and previews, and updates `Source/figures/workload_fct_comparison.pdf`. Missing Meta Hadoop baseline measurements remain absent.
