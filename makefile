all: DCC_cycloid.pdf clean view

DCC_cycloid.pdf: DCC_cycloid.tex DCC_cycloid.bib
	pdflatex DCC_cycloid.tex
	bibtex DCC_cycloid
	pdflatex DCC_cycloid.tex
	pdflatex DCC_cycloid.tex

view: DCC_cycloid.pdf
	open DCC_cycloid.pdf

clean:	
	rm -f DCC_cycloid.aux
	rm -f DCC_cycloid.log
	rm -f DCC_cycloid.out
	rm -f DCC_cycloid.bbl
	rm -f DCC_cycloid.blg
	rm -f DCC_cycloid.synctex.gz
	rm -f DCC_cycloid.fls
	rm -f DCC_cycloid.fdb_latexmk