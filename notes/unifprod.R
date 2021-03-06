library(Rcpp)

# see tests/test-unifprod.R for more details.

# Discretization of u(a,b) -> \int_a^b (u(a,t) u(t,b) - u(a,b)) dt
#  so output is y[i,k] = sum( x[i,i:(k-1)]*x[(i+1):k,k] ) - (k-1)*x[i,k]
#
# The key step is to do something like:
#     1 | x y z
#   a 2 | . . u
#     3 | . . v
#       +------
#     b   1 2 3
# --> [1,3]th entry of output is x*u + y*v - 2*z
# the first part is the inner product of a partial row with a partial column

# makes unifprod(), which works on matrices
unifprod.cpp <- paste(readLines("unifprod.cpp"),collapse="\n")

# makes unifprod_ut(), which works on vectors derived from column-ordered 
# upper-triangular portions of matrices, with diagonals
unifprod.ut.cpp <- paste(readLines("unifprod_ut.cpp"),collapse="\n")

# and does the two-argument version; see tests/test-unifprod.R.
biprod.ut.cpp <- paste(readLines("biprod_ut.cpp"),collapse="\n")


cppFunction(unifprod.cpp)
cppFunction(unifprod.ut.cpp)
cppFunction(biprod.ut.cpp)


