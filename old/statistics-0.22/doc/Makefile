figures = $(wildcard *.eps)

all: manual.pdf index.xhtml

manual.pdf: manual.texinfo
	tex manual.texinfo
	dvipdf manual

index.xhtml: manual.texinfo $(patsubst %.eps,%.png,$(figures))
	makeinfo --html manual.texinfo
	python html2xhtml.py index.html > index.xhtml
	rm index.html

$(patsubst %.eps,%.png,$(figures)): %.png: %.eps
	gs -r96 -dEPSCrop -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=png16m -sOutputFile=- -dBATCH -dNOPAUSE -q $< > $@


clean:
	rm -rf *.dvi *.aux *.log *.fn *.ky *.pg *.tp *.vr *.cp *.toc *.pdf *.txt *.png
