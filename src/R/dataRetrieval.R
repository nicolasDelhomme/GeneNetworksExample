#' ---
#' title: "Data retrieval"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(purrr)
  library(readr)
  library(tibble)
})

#' Data
#' First dataset:  #https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159151

dir.create(here("tmp"))

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE159151&format=file",
  destfile = here("tmp/GSE159151_RAW.tar")
)

files <- untar(here("tmp/GSE159151_RAW.tar"), exdir = here("tmp"))

files <- dir(here("tmp"), pattern = "*.txt.gz$", full.names = TRUE)

counts <- purrr::reduce(lapply(files, function(f){
  read_tsv(f,show_col_types = FALSE) %>% arrange(desc(2)) %>% 
    filter(!duplicated(`Gene symbol`))}),full_join)

#' Second dataset: 
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE177054&format=file&file=GSE177054%5FGene%5FExpression%2Exlsx",
  destfile = here("tmp/GSE177054.xlslx")
)

counts <- readxl::read_xlsx(here("tmp/GSE177054.xlslx")) %>% 
  mutate(avg=rowMeans(select(.,starts_with("count")))) %>% 
  arrange(desc(avg)) %>% filter(!duplicated(gene_name)) %>% 
  select("gene_name",starts_with("count")) %>% 
  full_join(counts,by=c("gene_name"="Gene symbol")) %>% 
  column_to_rownames("gene_name") %>% as.matrix()

#' # Export
dir.create(here("analysis/gene_network"),showWarnings=FALSE)
save(counts,file=here("analysis/gene_network/raw_counts.rda"))

#' # Cleanup
file.remove(dir(here("tmp"),full.names=TRUE))
file.remove(here("tmp"))

#' # Session Info
#' ```{r session info,echo=FALSE}
#' sessionInfo()
#' ```
