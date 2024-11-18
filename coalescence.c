#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "view.h"

#ifndef LEVEL
# define LEVEL 8
#endif

//scalar f[], * interfaces = {f};

double R1 = 0.4; //set the radius for the left droplet (R1<1.)
double R2 = 0.5; //set the radius for the right (R2<1.)
double muc1 = 1.8e-5; //set outer fluid dynamic viscosity
double muc2 = 1.e-3; //set inner fluid dynamic viscosity
double rhoc1 = 1.; //set outer fluid density
double rhoc2 = 1000.; //set inner fluid density
double sigmac = 5; //set surface tension coefficient

double uvelc = 1.; //set colliding speed (uniform velocity ??)
double runtime = 6.; //set runtime length
//#include "two-phase-generic.h"



int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  init_grid(1 << LEVEL); // Higher grid resolution
  const face vector muc[] = {muc1, muc2}; // {outer field, inner field}
  const face vector rhoc[] = {rhoc1, rhoc2}; // Density

  face vector rho[]; // Explicitly define rho as a face vector

  mu = muc;
  rho = rhoc;
  f.sigma = sigmac;
  run();
}

event init (t = 0)
{
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(R1)),
		    - (sq(x - 1.) + sq(y) - sq(R2))));
  double uvel = uvelc;
  foreach() {
    if (f[] > 1e-6){ //Apply velocity only within the droplets
      u.x[] = - sign(x)*f[] * uvel; //how to assign velocity to each droplet?
    }
   boundary({f,u});
  }
}

event movie (t += 0.04; t <= runtime)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");
  box();
  save ("movie.mp4");
}
