#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "contact-embed.h"
#include "view.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

int LEVEL = 7;

double theta_deg = 30.0; // gravity angle in degrees
double c_angle =   85.0; // contact angle in degrees
double gravity =   9.81; // gravity magnitude in m/s^2
double L =         1.0;  // domain size in meters
double t_end =     30.0; // simulation time in seconds
char output_dir[80];     // directory for output files
scalar f0;

int main(int argc, char * argv[]) {
  if (argc > 1)
    c_angle = atof(argv[1]);
  if (argc > 2)
    theta_deg = atof(argv[2]);
  if (argc > 3)
    gravity = atof(argv[3]);

  sprintf(output_dir, "output-g_angle=%.2g-c_angle=%.2g", theta_deg, c_angle);
  mkdir(output_dir, 0777);

  size (L);
  origin (-L/2., -L/2.);
  init_grid(1 << LEVEL);

  // Water below (denser, more viscous), air above
  rho1 = 1000.0;           // density of water [kg/m^3]
  rho2 = 1.225;            // density of air [kg/m^3]
  mu1 = 1e-3;              // dynamic viscosity of water [Pa.s] = [kg/(m·s)]
  mu2 = 1.8e-5;            // dynamic viscosity of air [Pa.s]

  f.sigma = 0.0728;        // surface tension water-air [N/m]

  const scalar c[] = c_angle * pi / 180.;
  contact_angle = c;

  run();
}

event acceleration (i++) {
  const double theta = theta_deg * pi/180.;
  face vector av = a;

  foreach_face(x)
    av.x[] = gravity * sin(theta);     // body force per unit mass [m/s^2]
  foreach_face(y)
    av.y[] = -gravity * cos(theta);    // body force per unit mass [m/s^2]
  a = av;
}

event init (t = 0) {
  // Initialize with lower half liquid, upper half gas,
  // Interface aligned perpendicular to gravity vector
  const double theta = theta_deg * pi/180.;
  //fraction (f,-y); 
  fraction(f,  +x*sin(theta) - y*cos(theta));
}


event stop (t = 20) {
  return 1;
  if (i == 0){
    f0 = f;
  }else{
    double max_diff = 0;

    foreach() {
      double df = fabs(f[] - f0[]);
      if (df > max_diff)
        max_diff = df;
    }
  
    if (max_diff < 1e-8)
      return 1;
  
    foreach()
      f0[] = f[];
  }
  
}

event output_interface (t+=0.5; last) {
  char filename[128];
  //sprintf(filename, "%s/interface.dat", output_dir);
  //FILE * fp = fopen(filename, "w");
  //output_facets(f, fp);
  //fclose(fp);

  sprintf(filename, "%s/movie.mp4", output_dir);
  view (fov = 20);
  draw_vof("f", lw = 2);
  squares("f", linear = true, min = 0, max = 1);

  char text[80];
  sprintf(text, "t = %.2f", t);
  draw_string(text, pos = 1, size = 20);
  save(filename);
}