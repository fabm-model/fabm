#
# Makefile for making the source code documentation.
#

SRCDIR = ../src

SHELL	   = /bin/sh

all: doc

doc:
	$(MAKE) -C $(SRCDIR) doc
	latex  report.tex
	bibtex report
	latex  report.tex
	latex  report.tex

report.dvi:
	$(MAKE) -C $(SRCDIR) doc
	latex  report.tex
	bibtex report
	latex  report.tex
	latex  report.tex

figures:
	wget -nd -nH -q http://www.gotm.net/data_files/gotm-figs-4.0.gz
	tar -xvzf gotm-figs-4.0.gz
	$(RM) gotm-figs-4.0.gz

clean:
	$(RM) *.log *.toc *.aux *.bbl *.blg *.bac *~

realclean: clean
	$(RM) *.dvi *.ps *.pdf
	$(RM) -rf html/
	$(RM) -rf $(PROTEXSRC)

distclean: realclean
	$(RM) *.gz *.zip *.Z
	$(RM) html html_scp pdf_scp
	$(RM) -rf figures

