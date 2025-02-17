## Library import
library(ggplot2)
library(reshape2)  # for melting the matrix
library(gridExtra)  # for arranging plots

library(vegan)
library(permute)
library(lattice)
library(dplyr)

# Set seed for reproducibility (optional)
set.seed(123)



### 1. Create binary matrix ###

# Generate a binary matrix
binary_matrix <- matrix(rbinom(20 * 10, 1, 0.5), nrow = 20, ncol = 10,
                        dimnames = list(paste0("agent_", LETTERS[1:20]),
                                        paste0("item_", 1:10)))
# Print the matrix
print(binary_matrix)



### 2. Order binary matrix ###

sort_matrix <- function(M) {
  # Sort columns by column sums in descending order
  M <- M[, order(colSums(M), decreasing = TRUE)]
  # Sort rows by row sums in ascending order
  M <- M[order(rowSums(M), decreasing = TRUE), ]
  return(M)
}

sorted_matrix <- sort_matrix(binary_matrix)
print(sorted_matrix)


### 3. Plot both matrix ###

# Function to reshape matrix for ggplot
plot_matrix <- function(mat, title) {
  df <- melt(mat)
  colnames(df) <- c("Agent", "Item", "Value")
  
  df$Agent <- as.numeric(factor(df$Agent))  # Convert Agent to numeric
  
  ggplot(df, aes(x = Item, y = Agent, fill = factor(Value))) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("0" = "cornsilk", "1" = "olivedrab"), 
                      name = "Value") +
    scale_y_reverse(breaks = unique(df$Agent), labels = levels(factor(df$Agent))) +  # Reverse y-axis
    labs(title = title, x = "Items", y = "Agents") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom")
}

# Plot original and sorted matrices
p1 <- plot_matrix(binary_matrix, "Original Matrix")
p2 <- plot_matrix(sorted_matrix, "Sorted Matrix (by row and column sums)")

# Arrange plots side by side
grid.arrange(p1, p2, ncol = 2)




### 4. Explore "Nestedtemp" functions from Vegan ###


# Documentation :"Function nestedchecker gives the number of checkerboard units, 
# or 2x2 submatrices where both species occur once but on different sites"
nestedchecker(binary_matrix)
nestedchecker(sorted_matrix)
# The matrix indeed doesn't need to be ordered.
# Must be a gaussienne x-axis = fill and middle at 0.5 fill


# Documentation : "Function nestedn0 implements nestedness measure N0 which is
# the number of absences from the sites which are richer than the most pauperate 
# site species occurs"
nestedn0(binary_matrix)
nestedn0(sorted_matrix)
# Result = 62, seems a lot



# Documentation : "Function nesteddisc implements discrepancy index which is the 
# number of ones that should be shifted to fill a row with ones in a table 
# arranged by species frequencies"

# Thought what was n0 (nb of holes in the 1s) -> different from the definition 
# from the oikos paper

# niter(number of iteration) is set to 200, let's leave it at that for now.
nesteddisc(binary_matrix)
nesteddisc(sorted_matrix)


## TEMP ##

# "Nestedtemp" automatically sorts the matrix, hence ordering it before 
# is not necessary.
nestedtemp(binary_matrix)
nestedtemp(sorted_matrix)
# When you run the whole code, the temp won't change but if you select only 
# these 2 lines, the temp change
# -> Mesure unstable, there must be randomness in the function 
# Gets stable with a seed 
# Doesn't do the same ordering (not the same result for binary and sorted matrix)


## NODF ##

# If set order = FALSE, the statistic is evaluated with the current matrix 
# Documentation :"ordering allowing tests of other meaningful hypothesis of 
# matrix structure than default ordering by row and column totals"

# First, the non ordered matrix but directly ordered in the function (order = TRUE)
nestednodf(binary_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)
# Binary + order = true same nodf than sorted + order = true


# "order" parameter is necessary when using an unsorted matrix 
# (take the non sorted matrix otherwise) :
nestednodf(binary_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)

# Already ordered matrix but ordered again in the function (order = TRUE)
nestednodf(sorted_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)

# Already ordered matrix but with order = FALSE 
nestednodf(sorted_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)
# Exactly the same result between these 2 functions with sorted_matrix



nestedbetasor(binary_matrix)
nestedbetasor(sorted_matrix)
# same here : no difference between sorted and unsorted matrix

nestedbetajac(binary_matrix)
nestedbetajac(sorted_matrix)
# same here : no difference between sorted and unsorted matrix

# Try with more iterations





### 5. Get temperature and plot temperature matrix ###

# Compute the nestedness nodf for the sorted matrix
nested_result_nodf <- nestednodf(sorted_matrix)
# S3 method for class "Nestednodf"
plot(nested_result_nodf, main = "Nestedness NODF (Sorted Matrix)")



# Plot the nestedtemp for sorted matrix with different colours
# Compute the nestedness temperature
nested_result_temp <- nestedtemp(sorted_matrix)

plot(nested_result_temp, kind = "temperature", col = rev(heat.colors(100)), 
     main = "Nestedness Temperature (Sorted Matrix)")
# The colours don't seem to make sense with the line
# Seems like the line we see on the plot is not exactly the same used 
# in the function



# Plot of nestedtemp for binary_matrix without different colours
out <- nestedtemp(binary_matrix)
out
plot(out, main = "Nestedness Temperature (Binary Matrix)")

plot(out, kind="incid", main = "Incidence Plot for Binary Matrix")
# Cross the line but surprising 





### 6. Temperature with item prevalence and agent inventory size ###

# Compute nestedtemp and transform it into a dataframe
out <- nestedtemp(binary_matrix)
my_df <- as.data.frame(out$u)
my_df$Agent <- rownames(my_df)

# Convert to long format
long_df <- reshape2::melt(my_df, id.vars = "Agent", variable.name = "Item", value.name = "Temperature")

# Compute metrics
item_prevalence <- colSums(binary_matrix)  # Item prevalence
agent_inventory <- rowSums(binary_matrix)  # Agent inventory size

# Add values to long_df
long_df$ItemPrevalence <- item_prevalence[as.character(long_df$Item)]
long_df$AgentInventory <- agent_inventory[as.character(long_df$Agent)]

# Find corresponding indices in binary_matrix
agent_indices <- match(long_df$Agent, rownames(binary_matrix))
item_indices <- match(long_df$Item, colnames(binary_matrix))

# Separate data based on binary_matrix values
long_df_1 <- long_df[binary_matrix[cbind(agent_indices, item_indices)] == 1, ]
long_df_0 <- long_df[binary_matrix[cbind(agent_indices, item_indices)] == 0, ]


# Plot for present items (1 = islands)
plot1 <- ggplot(long_df_1, aes(x = ItemPrevalence, 
                               y = AgentInventory, color = Temperature)) +
  geom_point(size = 3, alpha = 0.6, 
             position = position_jitter(width = 0.25, height = 0.25)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) + 
  labs(title = "Temperature - Islands Items",
       x = "Item Prevalence",
       y = "Agent Inventory Size",
       color = "Temperature") +
  theme_minimal()

# Plot for absent items (0 = holes)
plot0 <- ggplot(long_df_0, aes(x = ItemPrevalence, 
                               y = AgentInventory, color = Temperature)) +   
  geom_point(size = 3, alpha = 0.6, 
             position = position_jitter(width = 0.25, height = 0.25)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(title = "Temperature - Holes Items",
       x = "Item Prevalence",
       y = "Agent Inventory Size",
       color = "Temperature") +
  theme_minimal()

# Arrange plots side by side
grid.arrange(plot1, plot0, ncol = 2)

# Do 2 different plots
#plot1
#plot0




### 7. OECOSIMU ###

# Since we don't want to do these tests with the phoible data, 
# Here is the same code as the phoible one but with our random data

# Question : Really r1 and c1 that we want ? 
# Seems to be r0 and c0 in the vegan package (diff from r00)



# 7.1. NODF

# NODF - r1 (p-value = greater)
out <- oecosimu(sorted_matrix, nestednodf, "r1", alt = "greater")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations - NODF (r1)")

# NODF - c0 (p-value = greater)
out <- oecosimu(sorted_matrix, nestednodf, "c0", alt = "greater")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations - NODF (c0)")


# 7.2. TEMP

# Temp - r1 (p-value = less)
out <- oecosimu(sorted_matrix, nestedtemp, "r1", alt = "less")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations - TEMP (r1)")

# Temp - c0 (p-value = less)
out <- oecosimu(sorted_matrix, nestedtemp, "c0", alt = "less")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations - TEMP (c0)")



# 7.3. Inverse matrix (to "add" c1)


# Since c1 doesn't seem to be in the vegan package, we can try the idea to
# invert the rows and columns so that we can apply r1 to this "new" matrix


# Create the inverted matrix

# Matrix inversion
sorted_matrix_inverted <- t(sorted_matrix) 
print(sorted_matrix_inverted)


# Convert the inverted matrix to long format
df_inverted <- melt(sorted_matrix_inverted)
colnames(df_inverted) <- c("Item", "Agent", "Value")

# Ensure Item and Agent are factors and reverse the order of Item
df_inverted$Item <- factor(df_inverted$Item, levels = rev(unique(df_inverted$Item)))
df_inverted$Agent <- factor(df_inverted$Agent, levels = unique(df_inverted$Agent))

# Create the plot
ggplot(df_inverted, aes(x = Agent, y = Item, fill = factor(Value))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "cornsilk", "1" = "olivedrab"), 
                    name = "Value") +
  labs(title = "Visualization of the sorted inverted matrix (for c1)",
       x = "Items",
       y = "Agents") +
  theme_minimal()


# Temp - r1 with inverted matrix (p-value = less)
out <- oecosimu(sorted_matrix_inverted, nestedtemp, "r1", alt = "less")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations (inverted matrix) - TEMP (r1)")

# NODF - r1 with inverted matrix (p-value = greater)
out <- oecosimu(sorted_matrix_inverted, nestednodf, "r1", alt = "greater")
out
densityplot(permustats(out), as.table = TRUE, layout = c(1,4), 
            main = "Density of permutations (inverted matrix) - NODF (r1)")



# Since we do not have c1, we can't verify that we have the same results 
# between the normal matrix with c1 and the inverted matrix with r1
