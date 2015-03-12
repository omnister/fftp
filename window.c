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
   double a;
   double b;
   double c;
   double d;
   double e;
} WINLIST;

WINLIST windows[] = {
   // name          WINTYPE   n, a        b        c         d       e
   { "blackman",    BLACK,    1, 0.42659, 0.49656, 0.076849, 0.0000, 0.0},
   { "bnut",        BNUT,     2, 0.36358, 0.48918, 0.136599, 0.0106, 0.0}, 
   { "nuttall",     NUTT,     1, 0.35577, 0.48739, 0.144232, 0.0126, 0.0},
   { "flat",        FLAT,     1, 1.0,     1.93,    1.29,     0.388,  0.032}, 
   { "hamming",     HAMM,     3, 0.54,    0.46,    0.00,     0.0000, 0.0},
   { "hanning",     HANN,     3, 0.5,     0.50,    0.00,     0.0000, 0.0},
   { "kaiser",      KAISER,   1, 0.0,     0.0,     0.0,      0.0000, 0.0},
   { "rectangular", RECT,     1, 1.00,    0.00,    0.00,     0.0000, 0.0},
   { "none", 	    NONE,     1, 0.0,     0.0,     0.00,     0.0000, 0.0}
};

// forward references
double I0(double b);		// bessel function order 0
double bfromsl(double sl);	// compute b parameter from sidelobe level
double kaiser_window(int N, int L);
double kaiser_rbw(double sidelobe, int L);

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


// compute RBW 
int wintype_rbw(WINTYPE type, double *rbw) {
    int i;
    double rbw_table;
    double a,b,c,d,e;
    double n;
    double N=1024;
    int found;
    double scale;
    double theta;


    // get window parameters from table
    found = 0;
    for (i=0; (windows[i].wintype != NONE) && !found ; i++) {
       if (windows[i].wintype == type) {
	  a = windows[i].a;
	  b = windows[i].b;
	  c = windows[i].c;
	  d = windows[i].d;
	  e = windows[i].e;
	  found++;
       } 
    } 
    if (!found) {
	return(-1);
    }

    double area=0.0;
    double height=0.0;

    for (n = 0; n < N; n++) {

	if (type == KAISER) {
	    scale = kaiser_window(n,N);
	} else {
	    theta = (double)n*2.0*M_PI/(N-1.0);
	    scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta)+e*cos(4.0*theta));
	}
	
	area += scale*scale;
	height += scale;
    }

    *rbw = N*area/(height*height);

    // fprintf(stderr, "window rbw:%g computed:%g\n", N*area/(height*height));

    return(0);
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

    WINTYPE type=BLACK;
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

// compute window coefficients, normalize gain to unity and put them in win[]
// for future use...
int window_get(double *win, int N, WINTYPE type) {
    int n;
    double a, b, c, d, e;
    double gain;
    double scale;
    double theta;
    double area;
    double avg;

    int i;
    int found=0;
     
    // fprintf(stderr, "#opts.singlesided=%d null=%d\n", singlesided, nulldc);

    // get window parameters from table
    for (i=0; (windows[i].wintype != NONE) && !found ; i++) {
       if (windows[i].wintype == type) {
	  a = windows[i].a;
	  b = windows[i].b;
	  c = windows[i].c;
	  d = windows[i].d;
	  e = windows[i].e;
	  found++;
       } 
    } 
    if (!found) {
	return(-1);
    }

    area=0.0;
    avg=0.0;;
    for (n = 0; n < N; n++) {

	if (type == KAISER) {
	    scale = kaiser_window(n,N);
	} else {
	    theta = (double)n*2.0*M_PI/(N-1.0);
	    scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta)+e*cos(4.0*theta));
	}

	win[n] = scale;		// copy the coefficients
        area += scale;		// compute window weighting factor
    }

    gain =  area/((double)N);
    
    // fprintf(stderr,"#reav=%g,imav=%g,area=%g, g=%g\n", reavg, imavg, area, gain);

    for (n = 0; n < N; n++) {
	win[n] *= (1.0/gain);	// normalize gain
    }

    return(0);
}


// take a complex array in[] and window it with real coefficients in win[]
// optionally null out DC term or compensate DC for singlesided display
int window_do(COMPLEX *in, double *win, int N, int nulldc, int singlesided) {
    int n;
    int i;
    double scale;
    double theta;
    double area;
    double reavg;
    double imavg;
     
    reavg=0.0;
    imavg=0.0;
    for (n = 0; n < N; n++) {
	in[n].re = in[n].re*win[n];	// window the data
	in[n].im = in[n].im*win[n];

        reavg += in[n].re;	// compute integrated real dc value
        imavg += in[n].im;	// compute integrated imag dc value
    }

    if (nulldc) { 		// null out DC term completely
        for (n = 0; n < N; n++) {
	    in[n].re -= (win[n]*reavg/(double)N);
	    in[n].im -= (win[n]*imavg/(double)N);
        }
    } else {
	if (singlesided) { 	// subtract a windowed DC term to adjust DC level
            for (n = 0; n < N; n++) {
	        // makes dc term factor sqrt(2) lower
	        in[n].re -= (1.0+1.0/sqrt(2.0))*(win[n]*reavg/((double)N));
	        in[n].im -= (1.0+1.0/sqrt(2.0))*(win[n]*imavg/((double)N));
	    }
	}
    }

    return(0);
}

int window(COMPLEX *in, int N, WINTYPE type, int nulldc, int singlesided) {
    int n;
    double a, b, c, d, e;
    double gain;
    double scale;
    double theta;
    double area;
    double reavg;
    double imavg;

    int i;
    int found=0;
     
    // fprintf(stderr, "#opts.singlesided=%d null=%d\n", singlesided, nulldc);

    // get window parameters from table
    for (i=0; (windows[i].wintype != NONE) && !found ; i++) {
       if (windows[i].wintype == type) {
	  a = windows[i].a;
	  b = windows[i].b;
	  c = windows[i].c;
	  d = windows[i].d;
	  e = windows[i].e;
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
	    scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta)+e*cos(4.0*theta));
	}

	in[n].re = in[n].re*scale;	// window the data
	in[n].im = in[n].im*scale;

        area+=scale;		// compute window weighting factor
        reavg += in[n].re;	// compute integrated real dc value
        imavg += in[n].im;	// compute integrated imag dc value
    }

    // FIXME, should use lazy evaluation for calculation window gain. 
    // Needs only to be done once instead of in inner loop. 

    gain =  area/((double)N);
    
    // fprintf(stderr,"#reav=%g,imav=%g,area=%g, g=%g\n", reavg, imavg, area, gain);

    for (n = 0; n < N; n++) {
	in[n].re *= (1.0/gain);	// normalize gain
	in[n].im *= (1.0/gain);	// normalize gain

        // FIXME: need to know whether single-sided or double-sided
	// only correct DC term for single-sided spectrum

	if (nulldc) {
	    // subtract a windowed DC term to adjust DC level
	    if (type == KAISER) {
		scale = kaiser_window(n,N);
	    } else {
		theta = (double)n*2.0*M_PI/(N-1.0);
		scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta)+e*cos(4.0*theta));
	    }

	    // null out DC term completely
                 
	    in[n].re -= (1.0/(gain*gain))*scale*reavg/((double)N);
	    in[n].im -= (1.0/(gain*gain))*scale*imavg/((double)N);

            // fprintf(stderr,"%d %g %d\n", n, in[n].re, N);
            // fprintf(stderr,"#scale %g, reavg %g, gain %g\n", scale, reavg, gain);
            // fprintf(stderr,"%d %g\n", n, scale);
	} else {
	    if (singlesided) {
		// subtract a windowed DC term to adjust DC level
		if (type == KAISER) {
		    scale = kaiser_window(n,N);
		} else {
		    theta = (double)n*2.0*M_PI/(N-1.0);
		    scale = (a-b*cos(theta)+c*cos(2.0*theta)-d*cos(3.0*theta)+e*cos(4.0*theta));
		}
                 
		// makes dc term factor sqrt(2) lower
		in[n].re -= (1.0+1.0/sqrt(2.0))*((1.0/(gain*gain))*scale*reavg/((double)N));
		in[n].im -= (1.0+1.0/sqrt(2.0))*((1.0/(gain*gain))*scale*imavg/((double)N));
	    }
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
