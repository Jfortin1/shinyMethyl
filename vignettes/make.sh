R CMD sweave shinyMethyl.Rnw
R CMD pdflatex shinyMethyl.tex
bibtex shinyMethyl.aux
R CMD pdflatex shinyMethyl.tex
R CMD pdflatex shinyMethyl.tex
open shinyMethyl.pdf


SWEAVE_OPTIONS="eval=FALSE" R CMD Sweave --engine=knitr::knitr --pdf shinyMethyl.Rnw