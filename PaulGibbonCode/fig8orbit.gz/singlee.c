
#include <stdio.h>
#include <math.h>

/*  Declare field and particle variables */

float Ey, Bz, Ex;       /* field components */
float a0;          /* normalised laser amplitude */
float phase;       /* laser phase */

float px, py;      /* particle momenta */
float vx, vy;     /* particle velocity */
float gam;      /* relativistic factor */
float x, y;       /* particle position */
float q_over_m;   /* particle charge/mass ratio */

float trun;       /* run time (normalised laser periods) */
float delta_t;     /*  timestep */
float time;       /* time variable */  
int nstep;        /*  # integration steps */

/* intermediate variables for particle-pusher routine */

float px_minus, py_minus, px_plus, py_plus, px_dash;
float beta, gam_n, t, s;

const float pi = 3.1415927; 

double sin();
double sqrt();

FILE *orbit_file;     /* file for writing out particle orbit */            

int main() {

    a0 = 1.0;                  
    trun = 3*2*pi;             
    delta_t = 0.1; 
    nstep = trun/delta_t;       

    /*  Open a file for storing the particle orbit */

    orbit_file = fopen("orbit.dat", "w");


    x = 0.0;             /* initial position of electron */
    y = 0.0;
    px = 0.0;            /* initial momentum of electron */
    py = 0.0;
    q_over_m = -1.0;    /* charge:mass ratio */

    time = 0.0;          /* zero time counter */


    /* ---------------- */
    /*  Main loop  **   */
    /* ---------------- */

    while (time < trun) {

/*   update laser field */

      phase = time - x;
      Ey = a0*sin(phase);
      Ex = 0.0;
      Bz = a0*sin(phase);

/*   Buneman 'particle pusher' algorithm: */

      beta = 0.5*q_over_m * delta_t;     /*  constant absorbing charge, mass and dt */

      px_minus = px + beta*Ex;            /* 1/2-acceleration */
      py_minus = py + beta*Ey;

      gam_n = sqrt(1.0 + px_minus*px_minus + py_minus*py_minus);

      t = beta*Bz/gam_n;               /*  intermediate variables */
      s = 2*t/(1.0 + t*t);

      px_dash = px_minus + py_minus*t;   
      py_plus = py_minus - px_dash*s;      /* rotation about B */
      px_plus = px_dash  + py_plus*t;     

      px = px_plus + beta*Ex;            /*  1/2-acceleration */
      py = py_plus + beta*Ey;

      gam = sqrt(1.0 + px*px + py*py);    /*  new relativistic factor */

/*  compute velocities */

      vx = px/gam;
      vy = py/gam;

/*  update particle position */

      x = x + vx*delta_t;
      y = y + vy*delta_t;

/*   write out positions & momenta to file */
    
/*     write (30,'(4f10.5)') x,y,px,py */
/*      printf("time: %f, x: %f, y: %f\n", time, x, y); */
      fprintf(orbit_file, "%f  %f  %f  %f\n", x, y, px, py); 

/*   increment time variable */

      time = time + delta_t;

    }


/*   close files */

    fclose(orbit_file);

}
