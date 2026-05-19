#test the density for rotated copula
library(copula)
theta <- 10^{-1}
dim = 4
clayton_cop <- claytonCopula(param = theta, dim = dim)  # Define the Clayton copula
U <- matrix(c(0.9, 0.7, 0.9, 0.4), nrow = 1)               # Input as a matrix
density <- dCopula(U,clayton_cop,log = T)               # Compute the density
print(density)

density <- dCopula(matrix(c(0.9, 0.7, 0.1,0.4), nrow = 1),clayton_cop,log = T)  
print(density)

comp_dict <- generate_index(dim)
comp=3
ll.c(U, theta, comp=comp)
U^{1-comp_dict[comp,]}*(1-U)^comp_dict[comp,]



