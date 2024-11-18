#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"
#include "tension.h"
#include "view.h"


scalar f1[], f2[], f[], * interfaces = {f1, f2};
// Using f1, f2, and f for conditional merging

double R1 = 0.3; //set the radius for the left droplet (R1<1.)
double R2 = 0.3; //set the radius for the right (R2<1.)
double muc1 = 1.8e-5; //set outer fluid dynamic viscosity
double muc2 = 1.e-3; //set inner fluid dynamic viscosity
double rhoc1 = 1.; //set outer fluid density
double rhoc2 = 1000.; //set inner fluid density
double sigmac = 5.0; //set surface tension coefficient
double offset = 0.0;  // Offset for off-center collisions

double uvelc = 1.0; //set colliding speed (uniform velocity ??)
double runtime = 6.; //set runtime length
double index = 1.5; // 1 for bouncing, 0 for merging
//#include "two-phase-generic.h"



int main()
{
  size (4.);
  origin (-L0/2., -L0/2.);
  //init_grid(1 << LEVEL); // Higher grid resolution
  const face vector muc[] = {muc1, muc2}; // {outer field, inner field}
  const face vector rhoc[] = {rhoc1, rhoc2}; // Density

  face vector rho[]; // Explicitly define rho as a face vector

  mu = muc;
  rho = rhoc;
  f1.sigma = f2.sigma = sigmac;
  run();
}

event init (t = 0)
{
  fraction (f1, - (sq(x + 1.) + sq(y) - sq(R1)));
  fraction (f2, - (sq(x - 1.) + sq(y + offset) - sq(R2)));
  fraction (f, max (- (sq(x + 1.) + sq(y) - sq(R1)),
		    - (sq(x - 1.) + sq(y) - sq(R2))));
//  foreach(){
//    double uvel = uvelc;
//    if (f1[] > 1e-6 && f2[] > 1e-6){ //Apply velocity only within the droplets
//      u.x[] = (f1[] - f2[])*uvel;
//    }
//  }

// Assign initial velocity towards each other
  foreach() {
    if (f1[] > 1e-6) {
      u.x[] = uvelc;
    }
    if (f2[] > 1e-6) {
      u.x[] = -uvelc;
    }
  }
}

// Check for coalescence condition and merge droplets if met
event coalesce (i++)
{
  double distance = 0.1;  // Define a distance threshold for coalescence
  double min_distance = HUGE; // Initialize with a large value

  // Calculate the minimum distance between droplets
  foreach() {
    if (f1[] > 1e-6 && f2[] > 1e-6) {
      double dx = x, dy = y; // Placeholder for center positions
      double dist = sqrt(sq(dx) + sq(dy));
      if (dist < min_distance) min_distance = dist;
    }
  }

  // Conditional merging of droplets based on distance threshold
  if (min_distance < distance) {
    // Merge f1 and f2 into a single field f
    foreach() {
      f[] = f1[] + f2[];  // Merge f1 and f2 into f for coalescence
      f1[] = f2[] = 0;     // Clear f1 and f2 after merging
    }
    interfaces = {f};      // Switch to the single interface field
    index = 0.;
  }
}



event movie (t += 0.04; t <= 6.)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  if (index < 1.) {
    draw_vof ("f");
  }
  if (index > 1.) {
    draw_vof ("f1");
    draw_vof ("f2");
  }
  box();
  save ("movie.mp4");
}

// Track droplet breakup and fragmentation
//event snapshot (t += 0.1)
//{
//  char fname[99];
//  sprintf (fname, "snapshot-%.2f.ppm", t);
//  output_ppm (f1, file = fname, n = 512);
//}

// Stop simulation at end time
//event end (t = runtime)
//{
//  printf ("Simulation finished at t = %g\n", t);
//}
