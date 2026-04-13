# Import necessary libraries
library(tidyverse)
library(httr)
library(jsonlite)

genes_list <- read.csv("data/processed/common_rbps_with_prot_name.csv")

get_omim_disease <- function(omim_id, api_key) {
  
  if (is.na(omim_id)) return(NA)
  
  url <- paste0(
    "https://api.omim.org/api/entry?mimNumber=",
    omim_id,
    "&include=geneMap",
    "&apiKey=", api_key,
    "&format=json"
  )
  
  res <- httr::GET(url)
  
  if (httr::status_code(res) != 200) return(NA)
  
  json_data <- jsonlite::fromJSON(
    httr::content(res, "text", encoding = "UTF-8")
  )
  
  tryCatch({
    
    entry <- json_data$omim$entryList[[1]]$entry
    
    # Extract phenotype info
    if (!is.null(entry$geneMap$phenotypeMapList)) {
      
      phenotypes <- sapply(entry$geneMap$phenotypeMapList, function(x) {
        x$phenotypeMap$phenotype
      })
      
      return(paste(unique(phenotypes), collapse = "; "))
    }
    
    return(NA)
    
  }, error = function(e) {
    return(NA)
  })
}

# Get Primary associated diseases
api_key <- "YOUR_OMIM_API_KEY"

genes_list$disease <- sapply(genes_list$correct_omim_id, get_omim_disease, api_key)



