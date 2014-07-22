NAME = aa

system := $(shell uname -s)

ifeq ($(system), Darwin)
OPDF := open
else
OPDF := evince
#OPDF := acroread
endif

PS2PDF = ps2pdf -dEPSCrop

pq:
	pdflatex $(NAME)
	$(OPDF) $(NAME).pdf

all p pdf $(NAME).pdf: $(NAME).tex
	pdflatex $(NAME)
	bibtex $(NAME)
	pdflatex $(NAME)
	pdflatex $(NAME)

pdfo po: $(NAME).pdf
	$(OPDF) $(NAME).pdf

qp:
	pdflatex $(NAME)
	$(OPDF) $(NAME).pdf

clean:
	-/bin/rm $(NAME).dvi $(NAME).aux *~ $(NAME).log $(NAME).out $(NAME).bbl $(NAME).d $(NAME).bak $(NAME).blg 
	-/bin/rm $(epsd) $(pdfd)
	-/bin/rm $(spdfd)

tclean: clean
	-/bin/rm $(NAME).pdf

# source figures:
xfigs := $(wildcard xfig/*.fig)
epss := $(wildcard eps/*.eps)
spdf := $(wildcard spdf/*.pdf)

# derived figures:
epsd  := $(xfigs:xfig/%fig=eps/%eps)
pdfd  := $(xfigs:xfig/%fig=pdf/%pdf) $(epss:eps/%eps=pdf/%pdf)
spdfd := $(spdf:spdf/%pdf=pdf/%pdf)

figs: $(pdfd)
	echo $(pdfd)
	/bin/cp -p $(spdf) pdf

eps/%.eps: xfig/%.fig
	(cd `pwd`/xfig && figtex2eps $(<:xfig/%=%))
	mv $(<:%.fig=%.eps) eps

pdf/%.pdf: eps/%.eps
	$(PS2PDF) $< $@
