// We solve the incompressible, variable-density, Navier–Stokes equations \
with interfaces and surface tension. We can solve either the axisymmetric\
 or planar version. We can used standard or “reduced” gravity. We also \
 test levelset interface tracking and a momentum formulation.

#define CASE2 1
//#define AXIS 1
//#define MOMENTUM 1
#define CLSVOF 1
// solved - [ERROR] sb value is not changing with time ???
// clsvof has to be run with levelset on
#define LEVELSET 1
//#define REDUCED 1


#if AXIS
# include "axi.h" // fixme: does not run with -catch
#endif

#if MOMENTUM
# include "momentum.h"
#else
#include "navier-stokes/centered.h"

#if CLSVOF
# include "two-phase-clsvof.h"
#elif LEVELSET
# include "two-phase-levelset.h"
#else
# include "two-phase.h"
#endif
#endif

#if LEVELSET
# include "integral.h"
#else
# include "tension.h"
#endif
#if REDUCED
# include "reduced.h"
#endif

#ifndef LEVEL
# define LEVEL 8
#endif


//The boundary conditions are slip lateral walls (the default) and no-slip on \
the right and left walls.
// n-normal, t-tangential components of B.C.s
// ntr or xyz in 3D
#if MOMENTUM
q.t[right] = dirichlet(0);
q.t[left]  = dirichlet(0);
#else
u.t[right] = dirichlet(0);
u.t[left]  = dirichlet(0);
#endif

int main() {
  //The domain will span [0:2]×[0:0.5] and will be resolved with\
  256×64 grid points.
  size (2 [1]); //How to specify [0:2]*[0:0.5] and the resolution?
  DT = 1. [0,1];
  //[T, L, M] - dimension analysis in the square braket
  init_grid (1 << LEVEL);

  //Hysing et al. consider two cases (1 and 2), with the densities, dynamic \
  viscosities and surface tension of fluid 1 and 2 given below.
  rho1 = 1000.[0], mu1 = 10.;  // works also with rho1 = [-3,0,1]
  // property for water
#if CASE2
  rho2 = 1., mu2 = 0.1;
#else
  rho2 = 100., mu2 = 1.;
#endif

#if LEVELSET
  #if CASE2
  const scalar sigma[] = 1.96;
  #else
  const scalar sigma[] = 24.5;
  #endif
  d.sigmaf = sigma;
#else // !LEVELSET
  #if CASE2
  f.sigma = 1.96;
  #else
  f.sigma = 24.5;
  #endif
#endif // !LEVELSET

  // We reduce the tolerance on the Poisson and viscous solvers to\
  improve the accuracy.
  TOLERANCE = 1e-5 [*];
  // * is pointer ??
#if REDUCED
  G.x = -0.98;
  Z.x = 1.;
#endif
  run();
}

//event - to have variable timestep iterations, similar to for loop syntax
// t=0 (starting time), t<=5(the condition which must be verified for the event\
 to carry on), t+=1(the iteration operator)
//this one : start at t=0
event init (t = 0) {
  //The domain is a rectangle. We only simulate half the bubble.
  mask (y > 0.5 ? top : none);
  //The bubble is centered on (0.5,0) and has a radius of 0.25.
#if LEVELSET
  foreach()
    d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
#else
  fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
#endif
}

//We add the acceleration of gravity.
#if !REDUCED
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
}
#endif

//A utility function to check the convergence of the multigrid solvers.
// mgstats type is a structure (struct) containing several fields related to \
the performance or convergence of a multigrid solver (commonly used in \
numerical solvers for partial differential equations).

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    // mg.i(iteration counts?), mg.resa(current residual, should be +ve)
    // mg.reab (might be a previous residual value)
    // mg.nrelax (number of relaxation steps)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}
      // %d print an integer; %g print a floating-point value
      //This part is a ternary (conditional) operator. It evaluates the \
      condition mg.resb > 0, and based on whether this is true or false, \
      it performs one of two actions:
      //If mg.resb > 0 is true, it calculates exp(log(mg.resb / mg.resa) / mg.i).
      //If mg.resb > 0 is false, it returns 0. (a floating-point zero).

//We log the position of the center of mass of the bubble, its velocity\
and volume as well as convergence statistics for the multigrid solvers.
// i++ means executes at every timestep
event logfile (i++) {
  double xb = 0., vb = 0., sb = 0.;
  foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
    double dv = (1. - f[])*dv();
#if MOMENTUM
    vb += q.x[]*dv/rho(f[]);
#else
    vb += u.x[]*dv;
#endif
    xb += x*dv;
    sb += dv;
  }
  static double sb0 = 0.;
  if (i == 0) {
    printf ("t sb -1 xb vb dt perf.t perf.speed\n");
    sb0 = sb;
  }
  printf ("%g %g %g %g %g %g %g %g ",
	  t, (sb - sb0)/sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
//  printf ("%g %g %g %g %g %g %g %g %g ",
//  	 t, sb, sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
// !!! vb (rise velocity)
// !!! sb (relative volume difference)
#if !MOMENTUM
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
#endif
  putchar ('\n');
  fflush (stdout); // start a new line and ouput is displayed immediatly
}
//At t=3 we output the shape of the bubble
event interface (t = 3.) {
  output_facets (f, stderr);
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
}
#endif
