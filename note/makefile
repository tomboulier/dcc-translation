all: DCC_translation.pdf clean view

DCC_translation.pdf: DCC_translation.tex DCC_translation.bib
	pdflatex DCC_translation.tex
	bibtex DCC_translation
	pdflatex DCC_translation.tex
	pdflatex DCC_translation.tex

view: DCC_translation.pdf
	open DCC_translation.pdf

clean:	
	rm -f DCC_translation.aux
	rm -f DCC_translation.log
	rm -f DCC_translation.out
	rm -f DCC_translation.bbl
	rm -f DCC_translation.blg
	rm -f DCC_translation.synctex.gz
	rm -f DCC_translation.fls
	rm -f DCC_translation.fdb_latexmk