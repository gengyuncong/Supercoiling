NAME := UsersGuide
REFSNAME := references

default: pdf

all:	pdf refs

pdf:
	pdflatex $(NAME)
	bibtex $(NAME)
	pdflatex $(NAME)
	pdflatex $(NAME)
	rm -f `find . -name "*.aux"` $(NAME).bbl $(NAME).blg $(NAME).lof $(NAME).log $(NAME).lot $(NAME).toc
	pdfopen $(NAME).pdf

refs:
	bibtex-cleanup references.lib references.bib journals.dat

source:
	latex-demacro . $(NAME) UsersGuide.tex er-outline er-abbrev er-math
