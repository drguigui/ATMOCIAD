SRC = atmociad
LATEX = pdflatex
BIBTEX=bibtex
IDX=makeindex
ANNEXE = appendix/
#ORGA = organisation/
PART1 = ../photon/resultat
PART2 = ../electron/resultat
# PART3 = part3/




RERUN = '(There were undefined citations|Rerun to get (cross-references|the bars) right)'
UNDEFINED = '((Reference|Citation).*undefined)|(Label.*multiply defined)'

.PHONY:all

all: initial ${SRC}.pdf view
propre:
	@echo "Suppression des fichiers de compilation"
	@rm -rf *.log *.aux *.dvi *.toc *.lot *.lof *.out *.brf
	@rm -f ${ANNEXE}*.aux  ${ORGA}*.aux ${PART1}*.aux ${PART2}*.aux ${PART3}*.aux
	@echo "Suppression des fichiers bibtex"
	@rm -rf *.bbl *.blg *.bak
complet:
	@echo "Compilation complete"
	@if [ ! -e version.tex ]; then echo '  \\newcommand{\myRevision}{0} ' >> version.tex  ; fi
	@${LATEX} ${SRC}
	@if [ -e ${SRC}.bib ]; then ${BIBTEX} ${SRC} && echo "salut les gars" ; fi
	@${IDX} ${SRC}.nlo -s nomencl.ist -o ${SRC}.nls
	@${LATEX} ${SRC}
	@${LATEX} ${SRC}


svn: 
	@echo "Commit SVN"
	svn commit
	@echo '\\newcommand{\myRevision}{' > version.tex
	@svnversion -n >> version.tex
	@echo "}" >> version.tex
	


#commit: complet propre svn

commit: initial svn complet

initial: propre
	 @echo "Suppression des fichiers cibles"
	 @rm -rf *.pdf *.ps

%.pdf:%.tex
#	@echo '\\newcommand{\myRevision}{' > version.tex ;	
#	@svnversion -n >> version.tex;	
#	@echo "}" >> version.tex; 
	${LATEX} $<
	@if [ -e ${SRC}.bib ]; then ${BIBTEX} $* ; fi
	@if egrep -q ${RERUN} $*.log ; then echo "Passe Latex" && ${LATEX} $< ; fi
	@if egrep -q ${RERUN} $*.log ; then echo "Passe Latex" && ${LATEX} $< ; fi
	@${IDX} ${SRC}.nlo -s nomencl.ist -o ${SRC}.nls
	${LATEX} $<
	${LATEX} $<
	@echo "Citations ou ref indéfinies : "
	@egrep -i ${UNDEFINED} $*.log || echo "Aucune"


clean:propre

view:
	@evince ${SRC}.pdf &
