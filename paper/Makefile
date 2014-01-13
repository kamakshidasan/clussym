# set latexfile to the name of the main file without the .tex
latexfile = symmetryr1

# if no references, run "touch latexfile.bbl; touch referencefile.bib"
# or delete dependency in .dvi and tarball sections
referencefile = db

TEX = pdflatex

# support subdirectories
#VPATH = Figs

# reruns latex if needed.  to get rid of this capability, delete the
# three lines after the rule.  Delete .bbl dependency if not using
# BibTeX references.
# idea from http://ctan.unsw.edu.au/help/uk-tex-faq/Makefile
PNGS=$(wildcard figures/*.png)
FIGS=$(wildcard xfig/*.fig)
PDFS=$(wildcard figures/*.pdf)
FIGPDFS=$(FIGS:.fig=.pdf)

%pdf_t %.pdf : %.fig
	fig2dev -L pdftex $< $@; \
	fig2dev -L pdftex_t -p $*.pdf $*.fig $*.pdf_t

$(latexfile).pdf : $(PDFS) $(FIGPDFS) $(PNGS)  $(latexfile).tex $(latexfile).bbl
		while ($(TEX) $(latexfile) ; \
				grep -iq "Rerun to" $(latexfile).log ) do true ; \
					done


# keep .eps files in the same directory as the .fig files
#%.eps : %.fig
#		fig2dev -L eps $< > $(dir $< )$@
#$(latexfile).pdf : $(latexfile).ps
#		ps2pdf $(latexfile).ps $(latexfile).pdf




$(latexfile).bbl : db.bib
	bibtex $(latexfile)
	


#$(latexfile).ps : $(latexfile).dvi
#		dvips $(latexfile)

#ps : $(latexfile).ps 

# make can't know all the sourcefiles.  some file may not have
# sourcefiles, e.g. eps's taken from other documents. 
#$(latexfile).tar.gz : $(figures) $(latexfile).tex $(referencefile).bib
#		tar -czvf $(latexfile).tar.gz $^

#tarball: $(latexfile).tar.gz

