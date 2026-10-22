# DuoGran TPDS manuscript

LaTeX source for the DuoGran paper in IEEE Computer Society journal format.

## Build

Install a LaTeX distribution with IEEEtran (including its bibliography style), pdfLaTeX, BibTeX, and latexmk, then run:

```bash
make
```

The generated manuscript is `Source/prolet.pdf`. It is not tracked in Git.

The repository includes the main manuscript (`Source/prolet.tex`), bibliography (`Source/reference.bib`), build files, and all nine figure PDFs referenced by the manuscript. The figure PDFs are build inputs and are tracked.
