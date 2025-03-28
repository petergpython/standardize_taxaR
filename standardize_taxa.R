### 

library(dplyr)
library(tidyr)
library(httr)
library(jsonlite)

# read dataset 
df = read.csv("dataset.csv", header = TRUE )

##### correcting manually these species names as it occurs in a lot of accessions
df$species <- gsub('z.mays', 'mays', df$species)
df$species <- gsub('o.sativa', 'sativa', df$species)
######

############ create the input for the taxa standardization  
df$subTaxa2 <- ifelse(is.na(df$subTaxa), "", df$subTaxa)
df$full_taxa_input <- trimws(paste(df$genus, df$species, df$subTaxa2))

# Fill NA values with an empty string
df$spAuthor <- ifelse(is.na(df$spAuthor), "", df$spAuthor)

# Strip leading and trailing whitespace from the 'spAuthor' column
df$spAuthor <- trimws(df$spAuthor)

# Concatenate 'full_taxa_input' with 'spAuthor'
df$full_taxa_input <- paste(df$full_taxa_input, df$spAuthor)

# Strip leading and trailing whitespace from the 'full_taxa_input' column
df$full_taxa_input <- trimws(df$full_taxa_input)

####################
# function to query API of https://verifier.globalnames.org
query_taxa_resolver <- function(taxa, sources = c('196')){
  if (is.character(taxa)){
    taxa_format <- gsub(" ", "+", taxa)
    URL <- paste0('https://verifier.globalnames.org/api/v1/verifications/', taxa_format, 
                  '?data_sources=', paste(sources, collapse = "|"), 
                  '&all_matches=false&capitalize=true&species_group=false&fuzzy_uninomial=false&stats=false&main_taxon_threshold=0.8')
    tryCatch({
      r <- GET(URL)
      result <- content(r, "text", encoding = "UTF-8")
      result_dict <- fromJSON(result)
    }, error = function(e){
      return(NULL)
    })
  } else {
    return("undetermined")
  }
  return(result_dict)
}

# function to extract the best results from the query search 
extract_best_result <- function(list_res){
  final <- list()
  for (i in list_res){
    # added to handle the case one of the results of the query is NULL
    if (is.null(i)) {
      final <- append(final, list(c('null', 'no_match', 'no_match', 'no_match', 'no_match')))
    } else if (!("names" %in% names(i))) {
      final <- append(final, list(c('null', 'no_match', 'no_match', 'no_match', 'no_match')))
    } else {
      match_type <- i["names"][[1]]["matchType"]
      if (match_type != "NoMatch"){
        input_name <-  i["names"][[1]]$name
        matched_name <-i["names"][[1]]$bestResult$matchedName
        output_name <- i["names"][[1]]$bestResult$currentName
        status <-  i["names"][[1]]$bestResult$taxonomicStatus
        final <- append(final, list(c(input_name, matched_name, match_type, status, output_name)))
      } else {
        final <- append(final, list(c( i["names"][[1]]$name, 'no_match', 'no_match', 'no_match', 'no_match')))
      }
    }}
  return(final)
}

# taxa list to be standardised
taxa_list <- unique(trimws(na.omit(df$full_taxa_input)))

# loop trough taxa list and query the API
result_queries <- list()
counter <- 0
for (i in taxa_list){
  print(paste(round(counter / length(taxa_list) * 100, 2), "%", i))
  best_result <- query_taxa_resolver(i, c('196'))
  result_queries <- append(result_queries, list(best_result))
  counter <- counter + 1
}

# extract best result from the result of the queries
res <- extract_best_result(result_queries)

# create a taxa dictionary to be used to map the standardised taxa to the taxa in the dataset
taxa_standardized_df <- as.data.frame(do.call(rbind, res))
colnames(taxa_standardized_df) <- c('input_name', 'matched_name', 'match_type', 'status', 'output_name')
taxa_dictionary <- setNames(taxa_standardized_df$output_name, taxa_standardized_df$input_name)

# add standardized taxa to the combined dataset
df$taxa_standardized <- taxa_dictionary[df$full_taxa_input]

# save dataset as csv file
df_save <- apply(df,2,as.character)
write.csv(df_save, 'combined_df_standardized_taxa13_mar.csv', row.names = FALSE)

# save dataset 
taxa_standardized_df_save <- apply(taxa_standardized_df,2,as.character)

# save table summarizing how each taxa was standardized
write.csv(taxa_standardized_df_save, 'standardization_table.csv', row.names = FALSE)
