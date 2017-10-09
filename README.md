# libBeyn: A C++ Implementation of Beyn's Contour-Integration Method for Nonlinear Eigenproblems

<span style="font-variant:small caps">libbeyn</span> is a C++
library that implements the contour-integration method for
nonlinear eigenproblems described by W-J. Beyn in this paper:

+Wolf-J&uuml;rgen Beyn, ``An integral method for solving nonlinear 
 eigenvalue problems.'' *Linear Algebra and its Applications* **436** 
 3839 (May 2012).
+DOI: [https://doi.org/10.1016/j.laa.2011.03.030](https://doi.org/10.1016/j.laa.2011.03.030)
+ArXiV: [https://arxiv.org/1003.1580](https://arxiv.org/1003.1580)

More specifically, using a modified version of Beyn's notation,
we consider nonlinear eigenproblems of the form

**T**(*z*) **v** = 0

where **T** is a *D&times;D*-matrix-valued function of a complex
variable *z* and **v** is an *D*-dimensional complex vector.
(Beyn's notation uses lower-case *m* for the dimension *D* and
does not use boldface for matrices or vectors.)

## API

The `libBeyn` API is extremely simple; it basically just exports
a single routine, plus some associated support routines.

### Initialize data structure

```C++
 BeynSolver *CreateBeynSolver(int D, int L);
```

Allocate, initialize, and return a `BeynSolver` structure
for a problem of dimension *D* in which we expect no more
than *L* eigenvalues to reside within the integration we specify.

### Solve nonlinear eigenproblem (circular contour)

```C++

 int BeynSolve(BeynSolver *Solver,
               BeynFunction UserFunction, void *UserData,
               cdouble z0, double R, int N);
```

Compute all eigenvalues inside a circular contour,
centered at `z0` with radius `R,` using `N`-point trapezoidal-rule
quadrature for contour integrals.
(Here `cdouble` is short for `std::complex<double>`).

The parameters are as follows:

+ `BeynSolver *Solver`

    + Data structure created by a prior call to `CreateBeynSolver().`

+ `BeynFunction UserFunction`
+ `void *UserData

    + User-supplied function that performs linear algebra involving the **T** matrix (see below).

+  `cdouble z0`
+  `double R
+  `int N

    + Together these parameters specify the integration contour (a circle
centered at `z0` with radius `R`) and the quadrature strategy (`N`-point
rectangular rule).

The return value of `BeynSolve` is the number of eigenvalues
identified within the given contour. The actual eigenvalues
and eigenvectors are stored within the `Solver` structure (see below).

### Solve nonlinear eigenproblem (elliptical contour)

```C++
 int BeynSolve(BeynSolver *Solver,
               BeynFunction UserFunction, void *UserData,
               cdouble z0, double Rx, double Ry, int N);
```

A variant of the above routine in which the integration contour
is an ellipse centered at `z0` with horizontal and vertical
radii `Rx`, `Ry`. 

## User function

The prototype of the user-supplied function passed to `BeynSolve` is

```C++
  void UserFunction(cdouble z, void *UserData, HMatrix *VHat);
```

where `VHat` is an *M&times;L* complex-valued matrix. The user's function
should overwrite `VHat` with `T(z) \ VHat`, i.e. the function
should right-multiply the matrix `VHat` by the inverse of **T**(*z*).

### Eigenvalues and Eigenvectors

The eigenvalues and eigenvectors are stored in the `Lambda`
and `Eigenvectors` fields of `Solver.`
If the call to `BeynSolve` returns an integer `K`, then:

    +`Data->Lambda->GetEntry(k)` returns the `k`th eigenvalue

    +`Data->Eigenvectors->GetEntry(m,k)`
      returns the `m`th component of the `k`th eigenvector.

# Sample problem: Beyn example 4.11

The `libBeyn` distribution includes a code named `tBeyn411`
that uses `libBeyn` to solve the nonlinear eigenproblem
discussed in Example 4.11 of the Beyn paper referenced above,
which in turn borrowed the example from this reference:

+Kressner, D. "A block Newton method for nonlinear eigenvalue problems."
 *Numerische Mathematik* **114** 355 (2009)

+DOI: [https://doi.org/10.1007/s00211-009-0259-x](https://doi.org/10.1007/s00211-009-0259-x)

# Mode solver for <span style="font-variant:small caps">scuff-em</span>

Also packaged with `libBeyn` is a first stab at a mode solver for
[SCUFF-EM](http://homerreid.github.io/scuff-em-documentation).
This code, called `scuff-spectrum,` inputs:

    + a [SCUFF-EM geometry file](http://homerreid.github.io/scuff-em-documentation/reference),

    + a file containing a succession of `(&omega;0, R, N)` tuples.

For each `(&omega;0, R, N)` tuple in the file, `scuff-spectrum` executes
Beyn's method---for a circular contour of radius `R,` centered at `omega0,`
using `N`-point rectangular-rule quadrature---to compute all resonance 
frequencies of the SIE impedance matrix **M**(&omega;)---lying within 
the contor. (A resonance frequency &omega; is simply a frequency at
which **M(&omega;)** has a nontrivial nullspace, in which case the
usual SIE equation **M** **s**=**finc** can have nontrivial
surface-current solutions **s** even for vanishing incident 
field **finc**.)
