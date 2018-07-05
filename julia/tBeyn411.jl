#
# julia code that uses my Julia implementation of Beyn's method
# to solve the nonlinear eigenproblem of Example 4.11 in Beyn's paper
# 
# Homer Reid 10/2017
#

include("BeynSolve.jl");

##################################################
# user-supplied function passed to BeynSolve() 
##################################################
function Beyn411(z, VHat, MInvVHat, M)

  D   =  size(M,1);
  DE  =  2.0*D - 4.0*z/(6.0*D);  # diagonal entry
  ODE = -1.0*D -     z/(6.0*D);  # off-diagonal entry

  M*=0.0;
  M[1,1]=DE;  M[1,2]=ODE;        # top row
  for d=2:D-1                    # interior rows (2,D-1)
    M[d,d]=DE;
    M[d,d-1]=M[d,d+1]=ODE;
  end
  M[D,D-1] = ODE;                # bottom row
  M[D,D]   = 0.5*DE + z/(z-1.0);

@printf("Howdage %e+%e\n",real(z),imag(z))
  MInvVHat = M\VHat
  MInvVHat

end

##################################################
##################################################
##################################################
function SolveBeyn411(; z0=150.0+0.0*im, Rx=148.0, Ry=148.0, L=10, N=50)

  D = 400;
  M = im*zeros(D,D);

@printf("Howdage %e+%e\n",real(z0),imag(z0))
  (Lambda,V)=BeynSolve(z0, Rx, Ry, (z,V,MIV)->Beyn411(z,V,MIV,M), D; L=L, N=N)

  for n=1:length(Lambda)
    @printf("%2i: {%+.12e, %+.12e} \n",n,real(Lambda[n]),imag(Lambda[n]));
  end

end
