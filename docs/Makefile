userguide.pdf: userguide.tex refs.bib
	pdflatex userguide.tex
	bibtex userguide
	pdflatex userguide.tex
	pdflatex userguide.tex
	pandoc --toc -s userguide.tex -o README.md
	cp README.md ../README.md

clean:
	rm userguide.pdf README.md
