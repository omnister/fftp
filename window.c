#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fftlib.h"
#include "window.h"

typedef struct winlist {
   char *name;
   WINTYPE wintype;
   int n;	// number of chars name is unique in
   double rbw;
   double a;
   double b;
   double c;
   double d;
} WINLIST;

WINLIST windows[] = {
   // name          WINTYPE      n, rbw   a        b        c         d
   { "blackman",    BLACKMAN,    1, 1.73, 0.42659, 0.49656, 0.076849, 0.0000},
   { "nuttall",     NUTTALL,     1, 1.37, 0.35577, 0.48739, 0.144232, 0.0126},  // gain?
   { "bnut",        BNUT,        2, 1.37, 0.36358, 0.48918, 0.136599, 0.0106},  // gain?
   { "kaiser",      KAISER,      1, 1.00, 0.0,     0.0,     0.0,      0.0000}, 	// RBW?
   { "vonhanning",  HANNING,     1, 1.50, 0.5,     0.50,    0.00,     0.0000},
   { "hanning",     HANNING,     3, 1.50, 0.5,     0.50,    0.00,     0.0000},
   { "hamming",     HAMMING,     3, 1.37, 0.54,    0.46,    0.00,     0.0000},
   { "rectangular", RECTANGULAR, 1, 1.00, 1.00,    0.00,    0.00,     0.0000},
   { "none", 	    NONE,        1, 0.0,  0.0,     0.0,     0.00,     0.0000}
};

// forward references
double I0(double b);		// bessel function order 0
double bfromsl(double sl);	// compute b parameter from sidelobe level
double kaiser_window(int N, int L);

// cached private variables
static double b = 0.0;	
static double sidelobe = 80.0;

// set and remember sidelobe level for windows that
// require it (currently only kaiser);

int win_set_sidelobe(double sl) {
    extern double sidelobe;
    extern double b;

    sidelobe = sl;
    b = bfromsl(sl); // cache b value 
};

// set rbw based on wintype, return -1 on error, 0 on success

int wintype_rbw(WINTYPE wintype, double *rbw) {
    int i;
    for (i=0; windows[i].wintype != NONE; i++) {
       if (windows[i].wintype == wintype) {
	  *rbw = windows[i].rbw;
	  return(0);
       } 
    } 
    return(-1);
}

// set name based on wintype, return -1 on error, 0 on success

int wintypetoname(WINTYPE wintype, char **name) {
    int i;
    for (i=0; windows[i].wintype != NONE; i++) {
       if (windows[i].wintype == wintype) {
	  *name = windows[i].name;
	  return(0);
       } 
    } 
    return(-1);
}


int wintype(char *wname) {
    int i;
    int found=0;
    char *name;

    WINTYPE type=BLACKMAN;
    for (i=0; windows[i].wintype != NONE; i++) {
       if (strncasecmp(wname, windows[i].name, windows[i].n)==0) {
	  found++;
          type=windows[i].wintype;
       }
    }
    if (!found) {
       fprintf(stderr,"couldn't match \"%s\" window, using %s\n",
       		wname, name); 
    }
    // printf("%s %d\n", wname, type);
    return(type);
}

int window(COMPLEX *in, int N, WINTYPE type, int nulldc) {
    int n;
    double a, b, c, d;
    double gain;
    double scale;
    double theta;
    double area;
    double reavg;
    double imavg;

    int i;
    int found=0;

    // get window parameters from table
    for (i=0; (windows[i].wintype != NONE) && !found ; i++) {
       if (windows[i].wintype == type) {
	  a = windows[i].a;
	  b = windows[i].b;
	  c = windows[i].c;
	  d = windows[i].d;
	  found++;
       } 
    } 
    if (!found) {
	return(-1);
    }

    area=0.0;
    reavg=0.0;;
    imavg=0.0;
    for (n = 0; n < N; n++) {

	if (type == KAISER) {
	    scale = kaiser_window(n,N);
	} else {
	    theta = (double)n*2.0*M_PI/(N-1.0);
	    scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta));
	}

	in[n].re = in[n].re*scale;
	in[n].im = in[n].im*scale;
        area+=scale;		// compute window weighting factor
        if (nulldc) {
	   reavg += in[n].re;	// compute average real dc value
	   imavg += in[n].im;	// compute average imag dc value
	}
    }

    // FIXME, should use lazy evaluation for calculation window gain. 
    // Needs only to be done once instead of in inner loop. 

    gain =  area/((double)N);

    for (n = 0; n < N; n++) {
	in[n].re = in[n].re*(1.0/gain);	// normalize gain
	in[n].im = in[n].im*(1.0/gain);	// normalize gain

	if (nulldc) {
	    // subtract a windowed DC term to adjust DC level
	    if (type == KAISER) {
		scale = kaiser_window(n,N);
	    } else {
		theta = (double)n*2.0*M_PI/(N-1.0);
		scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta));
	    }
	    in[n].re -= ((1.0/gain)*scale*reavg/((double)N));
	    in[n].im -= ((1.0/gain)*scale*imavg/((double)N));

	    // makes dc term factor sqrt(2) lower
	    // in[n].re -= (1.0+1.0/sqrt(2.0))*((1.0/gain)*scale*reavg/((double)N));
	    // in[n].im -= (1.0+1.0/sqrt(2.0))*((1.0/gain)*scale*imavg/((double)N));
	} 
    }

    return(0);
}

void maintest() {
   wintype("hello");
   wintype("black");
   wintype("han");
   wintype("ha");
   wintype("ham");
   wintype("rect");
   wintype("kaiser");
}

double kaiser_window(int N, int L) {
  double a;
  extern double b;

  a=(L-1)/2.0;
  return (I0(b*sqrt(1.0-pow((((double)N-a)/a),2.0)))/I0(b));
}


// compute modified order zero bessel function
//
double I0(double b) {
   double n=0.0;
   double s=1.0;
   double u=1.0;

   do {
      n++;
      u*=pow(b/(2.0*n),2.0);
      s+=u;
   } while (fabs(u/s) > 1e-8);
   return (s);
}

// compute kaiser b value for given sidelobe level
// spinlab.wpi.edu/courses/ece503_2014/12-5kaiser_window_design.pdf
//
double bfromsl(double sl) {	
   if (sl > 120.0) sl=120.0;	// silently truncate sl

   if (sl <= 13.26) {
      return 0.0;
   } else if (sl > 13.26 && sl <= 60.0) { 
      return (0.76608*pow((sl-13.26),0.4)+0.09834*(sl-13.26));
   } else if (sl > 60.0 && sl <= 120.0) {
      return (0.12438*(sl+6.3));
   } 
}
