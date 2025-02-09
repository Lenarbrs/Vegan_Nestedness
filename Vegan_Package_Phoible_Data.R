# Vegan Package Description
# P.51 : null model (r1, c0, curveball, etc.)
# P. 145 : nested functions
# P. 153 : oecosimu (p value + simulations)

# Libraries
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)

library(vegan)
library(permute)
library(lattice)

# Upload Data
setwd("C:/Users/Léna/Desktop/Stage Urop/Data")
df_values <- read.csv("cldf-datasets-phoible-f36deac/cldf/values.csv", sep=",", header=TRUE)
df_languages <- read.csv("cldf-datasets-phoible-f36deac/cldf/languages.csv", sep=",", header=TRUE)




# 1. Select a random family (Alexey's code)
generate_binary_matrix <- function(languages, values) {
  n_langs <- table(languages$Family_Name)
  n_langs <- sort(n_langs[n_langs >= 4], decreasing = TRUE)
  
  if (length(n_langs) == 0) {
    stop("No family with at least 4 languages found.")
  }
  
  # Select a random family
  selected_family <- sample(names(n_langs), 1)
  
  # Select all languages from the selected family
  selected_langs <- languages$ID[languages$Family_Name == selected_family]
  
  # Select random contributions
  selected_contributions <- values %>%
    filter(Language_ID %in% selected_langs) %>%
    group_by(Language_ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Extract phonemes for selected contributions
  phonemes <- unique(values$Value[values$Contribution_ID %in% selected_contributions$Contribution_ID])
  phoneme_index <- setNames(seq_along(phonemes), phonemes)
  
  # Create Binary matrix 
  matrix_data <- matrix(0, nrow = length(selected_langs), ncol = length(phonemes),
                        dimnames = list(selected_langs, names(phoneme_index)))
  
  for (i in seq_along(selected_langs)) {
    lang_id <- selected_langs[i]
    phoneme_list <- values$Value[values$Contribution_ID == selected_contributions$Contribution_ID[selected_contributions$Language_ID == lang_id]]
    matrix_data[i, phoneme_list] <- 1
  }
  
  return(list(matrix = matrix_data, family = selected_family))
}

# 2. Sort the matrix by columns and rows sum
sort_matrix <- function(M) {
  if (ncol(M) > 1) {
    M <- M[, order(colSums(M), decreasing = TRUE)]
  }
  if (nrow(M) > 1) {
    M <- M[order(rowSums(M), decreasing = TRUE), ]
  }
  return(M)
}



# 3. Function to apply nestedness indices and display the results
apply_nestedness_metrics <- function(binary_matrix, sorted_matrix) {
  # Nestedtemp (for the sorted matrix)
  cat("\nNestedtemp results on the sorted matrix:\n")
  nested_result_temp_sorted <- nestedtemp(sorted_matrix)
  print(nested_result_temp_sorted)
  
  # Plot the temperature matrix
  plot(nested_result_temp_sorted, kind = "temperature", col = rev(heat.colors(100)), main = "Nested Temperature - Sorted Matrix")
  plot(nested_result_temp_sorted, kind="incid", main = "Incidence Matrix - Sorted Matrix")
  
  # Nestednodf (for the sorted matrix)
  cat("\nNestednodf results on the sorted matrix (order TRUE):\n")
  nested_result_nodf_sorted_true <- nestednodf(sorted_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)
  print(nested_result_nodf_sorted_true)
  
  cat("\nNestednodf results on the sorted matrix (order FALSE):\n")
  nested_result_nodf_sorted_false <- nestednodf(sorted_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)
  print(nested_result_nodf_sorted_false)
  
  # Nestedtemp for the binary matrix
  cat("\nNestedtemp results on the binary matrix:\n")
  nested_result_temp_bin <- nestedtemp(binary_matrix)
  print(nested_result_temp_bin)
  
  # Nestednodf for the binary matrix
  cat("\nNestednodf results on the binary matrix (order TRUE):\n")
  nested_result_nodf_bin_true <- nestednodf(binary_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)
  print(nested_result_nodf_bin_true)
  
  cat("\nNestednodf results on the binary matrix (order FALSE):\n")
  nested_result_nodf_bin_false <- nestednodf(binary_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)
  print(nested_result_nodf_bin_false)
}




# 4. Function to display matrices
plot_matrix <- function(mat, title) {
  df <- melt(mat)
  colnames(df) <- c("Language", "Phoneme", "Value")
  
  df$Language <- factor(df$Language, levels = rev(rownames(mat)))  # Keep the order of languages
  
  ggplot(df, aes(x = Phoneme, y = Language, fill = factor(Value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "cornsilk", "1" = "olivedrab"), name = "Presence") +
    labs(title = title, x = "Phonemes", y = "Languages") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom")
}





# 5. Generate and display the matrix
set.seed(123)
data <- generate_binary_matrix(df_languages, df_values)
binary_matrix <- data$matrix
sorted_matrix <- sort_matrix(binary_matrix)



# 6. Verify the matrices
cat("\nGenerated binary matrix:\n")
print(binary_matrix)

cat("\nGenerated sorted matrix:\n")
print(sorted_matrix)



# 7. Apply the nestedness indices
apply_nestedness_metrics(binary_matrix, sorted_matrix)

# Same result for nodf between: the already sorted matrix 
# (whether order = TRUE or FALSE in the function) 
# + the non-sorted matrix but with order = TRUE
# But lower nodf if the matrix is non-sorted and order = FALSE

# Different result between the already sorted matrix or not for temp:
# Sorted matrix: 38.24539 (Salishan family)
# Non-sorted matrix: 38.31987 (Salishan family)
# Not the same ordering between ours and the temp function's (and for some families seems better)



# 8. Plot the matrices
p1 <- plot_matrix(binary_matrix, paste("Original Matrix -", data$family))
p2 <- plot_matrix(sorted_matrix, paste("Sorted Matrix -", data$family))

# Display the matrices side by side
grid.arrange(p1, p2, ncol = 2)




# 9. OECOSIMU

# 9.1. NODF

# NODF - r1 (p-value = greater)
out <- oecosimu(sorted_matrix, nestednodf, "r1", alt = "greater")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), main = "Density of permutations - NODF (r1)")

# NODF - c0 (p-value = greater)
out <- oecosimu(sorted_matrix, nestednodf, "c0", alt = "greater")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), main = "Density of permutations - NODF (c0)")


# 9.2. TEMP

# Temp - r1 (p-value = less)
out <- oecosimu(sorted_matrix, nestedtemp, "r1", alt = "less")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), main = "Density of permutations - TEMP (r1)")

# Temp - c0 (p-value = less)
out <- oecosimu(sorted_matrix, nestedtemp, "c0", alt = "less")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), main = "Density of permutations - TEMP (c0)")

# c1 doesn't seem to exist, need to check



