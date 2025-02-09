## Library import
library(ggplot2)
library(reshape2)  # for melting the matrix
library(gridExtra)  # for arranging plots

library(vegan)
library(permute)
library(lattice)

# Set seed for reproducibility (optional)
set.seed(123)



### 1. Create binary matrix ###

# Generate an 8x20 binary matrix
binary_matrix <- matrix(rbinom(20 * 8, 1, 0.2), nrow = 20, ncol = 8,
                        dimnames = list(paste0("agent_", LETTERS[1:20]),
                                        paste0("item_", 1:8)))
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


# Function nestedchecker gives the number of checkerboard units, or 2x2 submatrices where both
# species occur once but on different sites
nestedchecker(binary_matrix)
nestedchecker(sorted_matrix)
# judging by the results, the matrix indeed doesn't need to be ordered.
# Must be a gaussienne x-axis = fill and middle at 0.5 fill

# Function nestedn0 implements nestedness measure N0 which is the number of absences from the
# sites which are richer than the most pauperate site species occurs 
nestedn0(binary_matrix)
nestedn0(sorted_matrix)
# Result = 62, seems a lot


# Function nesteddisc implements discrepancy index which is the number of ones that should be
# shifted to fill a row with ones in a table arranged by species frequencies 
# Thought what was n0 (nb of holes), different from the definition from the oikos paper
# niter(number of iteration) is set to 200, let's leave it at that for now.
nesteddisc(binary_matrix)
nesteddisc(sorted_matrix)

# "Nestedtemp" automatically sorts the matrix, hence ordering it before is not necessary.
nestedtemp(binary_matrix)
nestedtemp(sorted_matrix)
# When you run the whole code, the temp won't change but if you select only these 2 lines, the temp change
# Mesure unstable, there must be randomness in the function 
# Gets stable with a seed 
# Doesn't do the same ordering (not the same result for binary and sorted matrix)
# Ordering of the function seems better with better temp


# If set order = FALSE, the statistic is evaluated with the current matrix 
# ordering allowing tests of other meaningful hypothesis of matrix structure 
# than default ordering by row and column totals

# First, the non ordered matrix but directly ordered in the function (order = TRUE)
nestednodf(binary_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)
# Binary + order = true same nodf than sorted + order = true


# "order" parameter is necessary when using an unsorted matrix (take the non sorted matrix otherwise) :
nestednodf(binary_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)

# Already ordered matrix but ordered again in the function (order = TRUE)
nestednodf(sorted_matrix, order = TRUE, weighted = FALSE, wbinary = FALSE)

# Already ordered matrix but with order = FALSE 
# "order" parameter makes weird things when not used (on a sorted matrix) :
# Problem resolved : our sorted matrix was nested on the bottom left
# But the NODF function started to count the rows in the top left but there were no 1s
nestednodf(sorted_matrix, order = FALSE, weighted = FALSE, wbinary = FALSE)


nestedbetasor(binary_matrix)
nestedbetasor(sorted_matrix)
# same here : no difference between sorted and unsorted matrix

nestedbetajac(binary_matrix)
nestedbetajac(sorted_matrix)
# same here : no difference between sorted and unsorted matrix




### 5. Get temperature and plot temperature matrix ###

# Compute the nestedness nodf
nested_result_nodf <- nestednodf(sorted_matrix)
# S3 method for class "Nestednodf"
plot(nested_result_nodf, main = "Nestedness NODF (Sorted Matrix)")

# When we order ourselves in step 2 using our sorting function, instead of the one included in the functions
# When we create our matrix in part 1, Agent A is at the top
# But when we plot in part 3, the y-axis is reversed, meaning Agent A is at the bottom
# We then realize that the plot represents the inverse of what we are doing; it also shows our sorted matrix upside down
# In part 1, the nested part is at the bottom left
# In the plot, the nested part is at the top left, change the y-axis



# Compute the nestedness temperature
nested_result_temp <- nestedtemp(sorted_matrix)

# Plot the nested matrix using the correct method
plot(nested_result_temp, kind = "temperature", col = rev(heat.colors(100)), main = "Nestedness Temperature (Sorted Matrix)")
# The colours don't seem to make sense with the line
# Not exactly the same matrix as in the section 6 




### 6. Example from Vegan documentation ###

# Matrix temperature
out <- nestedtemp(binary_matrix)
out
plot(out, main = "Nestedness Temperature (Binary Matrix)")

plot(out, kind="incid", main = "Incidence Plot for Binary Matrix")
# Cross the line but surprising 

## Use oecosimu to assess the non-randomness of checker board units
nestedchecker(binary_matrix)
oecosimu(binary_matrix, nestedchecker, "quasiswap")
# Quasiswap
# + add alternative greater (test "one-sided")
## Another Null model and standardized checkerboard score
oecosimu(binary_matrix, nestedchecker, "r00", statistic = "C.score", alternative = "greater")
