library(here)

ptm <- proc.time()

rmarkdown::render(here('bone_marrow/analyze_bone_marrow.Rmd'))

rm(list = ls())

rmarkdown::render(here('cord_blood/analyze_cord_blood.Rmd'))

rm(list = ls())

proc.time() - ptm