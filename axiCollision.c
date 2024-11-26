#define AXIS 1
#define ADAPT 1
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

# include "tag.h"

#ifndef LEVEL
//# define LEVEL 8
# define LEVEL 10
#endif
#include "view.h"
//#include "pfplib.h"


//u.n[right] = neumann(0.);
//u.n[left] = neumann(0.);
u.t[left] = dirichlet(0.);
u.t[right] = dirichlet(0.);
//p[right] = dirichlet(0.);
//p[left] = dirichlet(0.);

//u.n[top] = neumann(0.);
//u.n[top] = dirichlet(0.);
//p[top] = dirichlet(0.);

int main(int argc, char *argv[]) {

size (10);
origin (-L0/2., 0.0);

DT = 1. [0,1];
init_grid (1 << LEVEL);

/*
  rho2 = 1000.[0], mu2 = 10.;  // works also with rho2 = [-3,0,1]
rho1 = 1., mu1 = 0.1;
*/ 
/* air and tetradecanedensity and viscosity ratios */
double rhoR, muR;
double We, Re, Oh;
We = 61.4;
Oh = 0.02643;

/* start pfplib for parameter readings */

// getPars(argc, argv);

/* dimensionless numbers 
if(dfetch("We", &We) != 1) {printf("missing We\n"); exit(1);}
if(dfetch("Oh", &Oh) != 1) {printf("missing Oh\n"); exit(1);}
	dfetch("Re", &Re);
*/

Re = sqrt(We)/Oh;

printf("We = %g, Re = %g, Oh = %g\n", We, Re, Oh);

/* rho and viscosity rations of two phases */
//rhoR=1.225/762.0; 
//muR=0.01827/2.128;
rhoR=1.225/762.0; 
//rhoR=.01; 
muR=0.01827/2.128;

rho1 = 1.[0]; mu1 = 1.0/Re;  // works also with rho2 = [-3,0,1]
rho2 = rho1*rhoR, mu2 = mu1*muR;
#if LEVELSET
  const scalar sigma[] = 5.0;
  d.sigmaf = sigma;
#else // !LEVELSET
  f.sigma = 1.0/We;
#endif // !LEVELSET

 TOLERANCE = 1e-6 [*];
  run();
}

event init (t = 0) {
//  mask (y > 5.0 ? top : none);
#if LEVELSET
  foreach()
    d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
#else
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(0.5)),
                    - (sq(x - 1.) + sq(y) - sq(0.5))));
#endif
  // initialize the velocity field
    double uvel = 0.5;
    foreach(){
    //  if(x > 1) u.x[] = -f[] * uvel;
    //  if(x <= 1) u.x[] =  f[] * uvel;
     u.x[] = -sign(x)* f[] * uvel;
    }
    
}
void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}
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
    //printf ("t sb -1 xb vb dt perf.t perf.speed\n");
    sb0 = sb;
  }
  if(pid() == 0) printf ("%g %g %g %g %g %g %g %g ",
	  t, (sb - sb0)/sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);
#if !MOMENTUM
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
#endif
  putchar ('\n');
  fflush (stdout);
}

event interface (t ={2., 2.1, 2.2, 2.3,2.4,2.5, 2.6, 2.7, 2.8, 2.9, 3.0}) {
//event interface (t ={0., .1}) {
//event interface (t ={0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22}) {
//  output_facets (f, stderr);
  remove_droplets(f, minsize = 20.0, bubbles = true);
}

event movie (t += 0.1; t <= 20.)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");
  box();
  save ("collision.mp4");
}

event snapshot (t = 1; t <= 20; t++)
{
  char name[80];
  sprintf (name, "dump-%03d", (int) t);
  dump (file = name);
}

#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){5e-6,1e-3,1e-3}, LEVEL);
}
#endif
