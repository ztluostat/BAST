### Redefine functions related to FEM ###
### code adapted from fdaPDE ###

# function to returns the index of the triangle containing the point
# returns multiple triangles for points on edges if all_tri == TRUE
# code adapted from fdaPDE
R_insideIndex2 <- function (mesh, location, all_tri = FALSE)
{
  #  insideIndex returns the index of the triangle containing the point
  # (X,Y) if such a triangle exists, and NaN otherwise.
  #  TRICOEF may have already been calculated for efficiency,
  #  but if the function is called with four arguments, it is calculated.
  
  
  eps=2.2204e-016
  small = 10000*eps
  
  nodes = mesh$nodes
  triangles = mesh$triangles
  X = location[1]
  Y = location[2]
  
  ntri   = dim(triangles)[[1]]
  indtri   = matrix(1:ntri,ncol=1)
  
  #  compute coefficients for computing barycentric coordinates if needed
  
  tricoef = fdaPDE:::R_tricoefCal(mesh)
  
  #  compute barycentric coordinates
  r3 = X - nodes[triangles[,3],1]
  s3 = Y - nodes[triangles[,3],2]
  lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
  lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
  lam3 = 1 - lam1 - lam2
  
  #  test these coordinates for a triple that are all between 0 and 1
  int  = (-small <= lam1 & lam1 <= 1+small) & 
    (-small <= lam2 & lam2 <= 1+small) & 
    (-small <= lam3 & lam3 <= 1+small)
  
  #  return the index of this triple, or NaN if it doesn't exist
  indi = indtri[int]
  if (length(indi)<1)
  {
    ind = NA
  }else if(all_tri) {
    ind = indi
  } else{
    ind = min(indi)
  }
  
  ind
}

