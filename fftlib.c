#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fftlib.h"


int snap2(int n) {	// snap to nearest power of two
    return (int)pow(2.0,(1.0+(int)(log((double)n-1.0)/log(2))));
}

COMPLEX *new_complex(int size) {
   return((COMPLEX *)malloc(sizeof(COMPLEX)*size));
}

// 1D Fast Fourier Transform (FFT).
// The array length must be a power of two
// from http://lux.vu/applets/frequency.pde

COMPLEX *fft_1d( COMPLEX *array, int n)
{
    double u_r, u_i, w_r, w_i, t_r, t_i;
    int ln, nv2, k, l, le, le1, j, ip, i;

    ln = (int) (log((double) n) / log(2.0) + 0.5);
    nv2 = n / 2;
    j = 1;
    for (i = 1; i < n; i++) {
	if (i < j) {
	    t_r = array[i-1].re;
	    t_i = array[i-1].im;
	    array[i-1].re = array[j-1].re;
	    array[i-1].im = array[j-1].im;
	    array[j-1].re = t_r;
	    array[j-1].im = t_i;
	}
	k = nv2;
	while (k < j) {
	    j = j-k;
	    k = k/2;
	}
	j = j + k;
    }

    for (l = 1; l <= ln; l++) {	/* loops thru stages */
	le = (int) (exp((double) l * log(2.0)) + 0.5);
	le1 = le / 2;
	u_r = 1.0;
	u_i = 0.0;
	w_r = cos(M_PI / (double) le1);
	w_i = -sin(M_PI / (double) le1);

	/* loops thru 1/2 twiddle values per stage */
	for (j = 1; j <= le1; j++) {	
	    /* loops thru points per 1/2 twiddle */
	    for (i = j; i <= n; i += le) {	
		ip = i + le1;
		t_r = array[ip-1].re * u_r - u_i * array[ip-1].im;
		t_i = array[ip-1].im * u_r + u_i * array[ip-1].re;

		array[ip-1].re = array[i-1].re - t_r;
		array[ip-1].im = array[i-1].im - t_i;

		array[i-1].re = array[i-1].re + t_r;
		array[i-1].im = array[i-1].im + t_i;
	    }
	    t_r = u_r * w_r - w_i * u_i;
	    u_i = w_r * u_i + w_i * u_r;
	    u_r = t_r;
	}
    }
    return array;
}				/* end of FFT_1d */

// Inverse FFT routine.
// The array length must be a power of two.

COMPLEX *ifft_1d(COMPLEX *array, int n)
{
    double u_r, u_i, w_r, w_i, t_r, t_i;
    int ln, nv2, k, l, le, le1, j, ip, i;

    ln = (int) (log((double) n) / log(2.0) + 0.5);
    nv2 = n / 2;
    j = 1;
    for (i = 1; i < n; i++) {
	if (i < j) {
	    t_r = array[i-1].re;
	    t_i = array[i-1].im;
	    array[i-1].re = array[j-1].re;
	    array[i-1].im = array[j-1].im;
	    array[j-1].re = t_r;
	    array[j-1].im = t_i;
	}
	k = nv2;
	while (k < j) {
	    j = j-k;
	    k = k / 2;
	}
	j = j + k;
    }

    for (l = 1; l <= ln; l++) {	/* loops thru stages */
	le = (int) (exp((double) l * log(2)) + 0.5);
	le1 = le / 2;
	u_r = 1.0;
	u_i = 0.0;
	w_r = cos(M_PI / (double) le1);
	w_i = sin(M_PI / (double) le1);
	/* loops thru 1/2 twiddle values per stage */
	for (j = 1; j <= le1; j++) {	
	    /* loops thru points per 1/2 twiddle */
	    for (i = j; i <= n; i += le) {	
		ip = i + le1;
		t_r = array[ip-1].re * u_r - u_i * array[ip-1].im;
		t_i = array[ip-1].im * u_r + u_i * array[ip-1].re;

		array[ip-1].re = array[i-1].re - t_r;
		array[ip-1].im = array[i-1].im- t_i;

		array[i-1].re = array[i-1].re + t_r;
		array[i-1].im = array[i-1].im + t_i;
	    }
	    t_r = u_r * w_r - w_i * u_i;
	    u_i = w_r * u_i + w_i * u_r;
	    u_r = t_r;
	}
    }
    return array;
}				/* end of ifft_1d */


/************Chirp Z Transform *****************
* N = # input samples. M = # output samples.
*
* fsam   = the sample frequency in Hz.
* fstart = the start frequency in Hz.
* fstop  = the stop frequency in Hz for the band over
*  which the transform is computed.
*
* (fstart-fstop)/M = new resolution  
*
* Note: this method returns an array of length L. 
* L = the returned transform length. L will always be
* larger than M.  The first M samples are the result.
*
* It is the callers responsibility to compute L, create
* g[L] and free the array after using it.
*
**********************************************/

int czt_blocksize(int N, int M) {
  return snap2(M+N);
}

COMPLEX *czt( 
	COMPLEX array[],	// input sequence, size N
        int N, 			// size of input sequence
	int M, 			// size of output sequence
	COMPLEX g[],		// array used to return results
	int L,			// size of scratch array
	double fstart, 		// start frequency
	double fstop,		// stop frequency
        double fsam)		// sampling frequency
{
    COMPLEX *h = new_complex(L);
    double theta0;
    double phi0;
    double psi;
    double a;
    double b;
    int n;
    int k;


    // fprintf(stderr,"n=%d m=%d fb=%g fe=%g fsam=%g\n", N,M,fstart,fstop,fsam);

    phi0 = 2.0 * M_PI * (fstop-fstart) / fsam / (M-1);
    theta0 = 2.0 * M_PI * fstart / fsam;

    /*** Create arc coefficients ***/
    for (n = 0; n < M; n++) {
	h[n].re = cos(n * n / 2.0 * phi0);
	h[n].im = sin(n * n / 2.0 * phi0);
    }
    for (n = M; n < L - N; n++) {
	h[n].re = 0.0;
	h[n].im = 0.0;
    }
    for (n = L - N; n < L; n++) {
	h[n].re = cos((L-n) * (L-n) / 2.0 * phi0);
	h[n].im = sin((L-n) * (L-n) / 2.0 * phi0);
    }

    /*** Prepare signal ***/
    for (n = 0; n < N; n++) {
	g[n].re = array[n].re;
	g[n].im = array[n].im;
    }
    for (n = N; n < L; n++) {
	g[n].re = 0.;
	g[n].im = 0.;
    }

    double s;
    double c;

    for (n = 0; n < N; n++) {
	psi = n * theta0 + n * n / 2.0 * phi0;
	c = cos(psi);
	s = -sin(psi);
	a = c*g[n].re - s*g[n].im;
	b = s*g[n].re + c*g[n].im;
	g[n].re = a;
	g[n].im = b;
    }
    g = fft_1d(g,L);		//fft of samples    
    h = fft_1d(h,L);		//fft of arc coeff

    // convolution in the time domain is 
    // multiplication in the frequency domain

    // multiplication in the time domain is convolution 
    // in the frequency domain 

    for (n = 0; n < L; n++) {
	c = g[n].re;
	s = g[n].im;
	a = c*h[n].re - s*h[n].im;
	b = s*h[n].re + c*h[n].im;

	g[n].re = a / L;	// for scaling purposes since
	g[n].im = b / L; 	// fft_1d does not use scale 
    }  

    free(h);

    g = ifft_1d(g,L);

    for (k = 0; k < M; k++) {
	psi = k*k/2.0 * phi0;
	c = cos(psi);
	s = -sin(psi);
	a = c*g[k].re - s*g[k].im;
	b = s*g[k].re + c*g[k].im;
	g[k].re = a;
	g[k].im = b;
    }

    return g;
}
