.PHONY: all clean

TEXBIN := $(HOME)/.local/share/texlive/2026/bin/x86_64-linux
export PATH := $(TEXBIN):$(PATH)
LATEXMK ?= latexmk

all:
	$(LATEXMK) -pdf -synctex=1 -interaction=nonstopmode -halt-on-error prolet.tex

clean:
	$(LATEXMK) -c prolet.tex
