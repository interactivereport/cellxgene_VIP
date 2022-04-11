#!/usr/bin/env Rscript

bookdown::render_book("index.Rmd", "bookdown::pdf_book")

# vi cellxgene_VIP.tex
# :g/includegraphics{https:\/\/interactivereport.github.io\/cellxgene_VIP\/tutorial\//s//includegraphics{/g

# Google and download missing latex style files to the same directory as .tex file or into this directory structure
#    /Users/bzhang1/Library/texmf/tex/latex
#    US-M095355:latex bzhang1$ tree
#    .
#    ├── algorithmic
#    │   └── algorithmic.sty
#    ├── authblk
#    │   └── authblk.sty
#    ├── chemgreek
#    │   └── chemgreek.sty
#    ├── framed
#    │   └── framed.sty
#    ├── fvextra
#    │   └── fvextra.sty
#    ├── mhchem
#    │   └── mhchem.sty
#    ├── multirow
#    │   └── multirow.sty
#    ├── subfigure
#    │   └── subfigure.sty
#    ├── threeparttable
#    │   └── threeparttable.sty
#    └── todonotes
#        └── todonotes.sty
#    


# Then run
# pdflatex -jobname=cellxgene_VIP cellxgene_VIP.tex 
