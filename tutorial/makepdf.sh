#!/usr/bin/env Rscript

bookdown::render_book("index.Rmd", "bookdown::pdf_book")

# vi cellxgene_VIP.tex
# :g/includegraphics{https:\/\/interactivereport.github.io\/cellxgene_VIP\/tutorial\//s//includegraphics{/g

# Google and download missing latex style files into this directory structure
#    /Users/bzhang1/Library/texmf/tex/latex
#    US-M095355:latex bzhang1$ tree




# Then run
# pdflatex -jobname=cellxgene_VIP cellxgene_VIP.tex 
