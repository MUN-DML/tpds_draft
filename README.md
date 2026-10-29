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

The script reads pinned records from the independent local `prolet-experiment-logs` archive. Set `PROLET_LOG_ARCHIVE` to its root on another machine. Raw FCT files are not included in this manuscript repository. It verifies input and FCT hashes, 128-host/80%/100-ms configurations, identical baseline flow identities and ideal FCTs, and complete per-scheme size-bin coverage. It regenerates the statistics snapshot, previews, and `Source/figures/workload_fct_comparison.pdf`. Both workloads include complete ConWeave and CONGA measurements. The original DuoGran curves are retained from archived raw logs whose SHA-256 hashes match the original plotted inputs exactly. Baseline curves use paired 100-ms inputs, whereas DuoGran retains its original 10-ms input and independent per-bin cohort; each duration and count is explicit in `plot_data.json`. Both AliStorage and Meta Hadoop use a linear lower range and a progressively logarithmic upper range with matching slopes at the transition. Workload-specific scales emphasize low and middle ranges to separate the curves. All ConWeave curves are fully visible with at least 5% vertical headroom for endpoint markers, and DuoGran is fully visible. CONGA points exceeding the displayed limits are clipped and recorded in `plot_data.json`; raw statistics remain unchanged.
