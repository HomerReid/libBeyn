##################################################
# julia code implementing Beyn's contour-integral
# method for nonlinear eigenproblems
# 
# homer reid 10/2017
##################################################

##################################################
# Implement Beyn's method to compute all eigenvalues
# of a DxD matrix lying within an elliptical contour
# centered at z0, with horizontal/vertical radius Rx/Ry,
# using N rectangular-rule quadrature points, assuming 
# the contour contains no more than L eigenvalues.
#                   
# The matrix is only ever referenced abstractly through
# the user-supplied UserBeynFunction, which has prototype
#                   
#  UserBeynFunction(z, VHat, MInvVHat)
#                   
# and which should set MInvVHat = M(z) \ VHat and 
# return MInvVHat.
#
# BeynSolve returns a tuple (Lambda,V)
# where Lambda[k] = kth eigenvalue,
#       V[:,k]    = kth eigenvector
##################################################
function BeynSolve(z0, Rx, Ry, UserBeynFunction, D; L=10, N=50)
   
  MInvVHat   = im*zeros(D,L);
  A0         = im*zeros(D,L);
  A1         = im*zeros(D,L);
  VHat       = randn(D,L) + im*randn(D,L);

  #################################################################
  # evaluate contour integrals
  #################################################################
  DeltaTheta = 2.0*pi/N;
  for n=0:N-1
    Theta      = n*DeltaTheta;
    CT         = cos(Theta);
    ST         = sin(Theta);
    z1         = Rx*CT + im*Ry*ST;
    dz         = (im*Rx*ST + Ry*CT)/N;
    MInvVHat   = UserBeynFunction(z0+z1, VHat, MInvVHat);
    A0        += dz*MInvVHat;
    A1        += z1*dz*MInvVHat;
  end
  
  #################################################################
  # linear-algebra postprocessing to extract eigenvalues/vectors 
  #################################################################
  (V,Sigma,W)=svd(A0);

  SigmaThreshold=1.0e-8;
  K=length(find( abs(Sigma) .> SigmaThreshold ) )
  @printf("Found %i nonzero singular values.\n",K);

  if (K==0)
    Lambda=V=im*zeros(0,0);
    return (Lambda,V)
  end

  if (K==L)
    @printf("** Warning: K=L=%i in BeynMethod (repeat with higher L)",K);
  end

  Vk = V[:,1:K];
  Sk = Sigma[1:K];
  Wk = W[:,1:K];
  B  = Vk' * A1 * Wk * diagm( 1.0 ./ Sk );

  (Lambda,V)=eig(B)

  Lambda += z0*ones(size(Lambda));

  Lambda,V

end

##################################################
# alternative entry point to BeynSolve for circular contour
# of radius R centered at z0 
##################################################
function BeynSolve(z0, R, UserBeynFunction, D; L=10, N=50)
  return BeynSolve(z0, R, R, UserBeynFunction, D; L=L, N=N)
end
