userguide.pdf: userguide.tex refs.bib
	pdflatex userguide.tex
	bibtex userguide
	pdflatex userguide.tex
	pdflatex userguide.tex
	pandoc --toc -t markdown -s userguide.tex -o README.md1
	pandoc README.md1 -t markdown-citations -o README.md --bibliography refs.bib
	cp README.md ../README.md

clean:
	rm userguide.pdf README.md
