#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>	// getopt
#include <getopt.h>

#include "fftlib.h"
#include "window.h"

#define MAXBUF 256	// for fgets()

#define OPTSN_DEFAULT	1024	// default block size

double dbparse(char *dbrefstring);
int parse_fspec(char *line, double *fstart, double *fstop);
void dumparray(COMPLEX *in, int npts, char *str);

typedef struct options {
    int n;
    double p;
    int wtype;
    double wbw;
    double rbw;
    double kaiser;
    int m;
    int fflag;
    double fstart;
    double fstop;
    int czt;
    int animate;
    int pdplot;
    int nodc;
    int singlesided;
    int peakhold;
    double samplerate;
    int dbmode;
    double dbref;
    int verbose;
    int in;
    int out;
    int spect;	// spectral density mode
} OPTS;

OPTS opts;

void init_opts() {
   opts.n = 0;
   opts.p = 1.0;
   opts.spect = 0;
   opts.wtype = 0;
   opts.wbw = 1.0;
   opts.rbw = 0.0;
   opts.kaiser = 120.0;		// kaiser sidelobe depth
   opts.m = 0;                  // CZT options
   opts.fflag = 0;		// 1=>fstart set, 2, fstop set, 0=> none
   opts.fstart = 0.0;           // CZT options
   opts.fstop = 0.0;            // CZT options
   opts.czt=0;                  // default to no CZT
   opts.animate = 0;
   opts.pdplot = 0;
   opts.nodc = 0;               // null dc (!implemented)
   opts.animate = 0;
   opts.singlesided = 0;
   opts.peakhold = 0;
   opts.samplerate = 1.0;
   opts.dbmode = 0;             // -db
   opts.dbref = 1.0;            // -db(m|v|u|c|<n.n>)
   opts.verbose=0;
   opts.in=0;
   opts.out=0;
}

void massage_opts() {
    int debug=1;
    char *name;

    if (wintype_rbw(opts.wtype, &opts.wbw) != 0) {
       fprintf(stderr, "couldn't get window bandwidth\n");
       exit(1);
    }

    if (wintypetoname(opts.wtype,  &name) != 0) {
       fprintf(stderr, "couldn't get window name\n");
       exit(1);
    }

    // chicken and egg problem:
    // if the user specifies the resolution bw with -b option
    // then we need the sample rate to compute N
    // How do we know sample rate until we read in the data?
    // How can we read in until we know N?
    // Solution: create a wrapper around fgets();
    // Then we can set our wrapper to save lines
    // in a queue when in pushback mode.  We can then 
    // read a small number (say 4) lines into a buffer
    // and returns the sample rate.  
    // This allows computation of N.  When
    // readin(N,buf) gets called later, for normal processing
    // we turnoff pushback mode. The existence
    // of the saved lines is checked and these data
    // points are read in as a special case and preread
    // data is cleared.

    // if opts.rbw !=0.0 then compute opts.n
    if (opts.rbw != 0.0) {

	COMPLEX readbuf[16];
	double sampletime;
	int ccflag=0;

	xfgets_pushback(1);	// turn on pushback
	readin(4,readbuf,&sampletime,0,&ccflag);
	xfgets_pushback(0); // turn off pushback
	// xfgets_dump();

	// fprintf(stderr, "sampletime is %g\n", sampletime);
	opts.samplerate = 1.0/sampletime;

        opts.n = (int)(opts.wbw*opts.samplerate)/opts.rbw;
    }

    if (opts.m != 0 && opts.n == 0) {
       opts.n = opts.m;
    } else {
       if (opts.n == 0) {
	   opts.n = OPTSN_DEFAULT;
       }
    }

    if(!(opts.czt)) {
        opts.in = snap2((int) (((double)opts.n)*opts.p));
        opts.out = opts.in;
	opts.m = opts.in;
    } else {
        if (opts.m == 0) {
           opts.m = (int) opts.p*opts.n;
        }
	// FOO
        opts.in = opts.n*opts.p;
        // opts.in = snap2(opts.n*opts.p);
        opts.out = snap2(opts.n*opts.p+opts.m);
    }

    if(opts.verbose) {
        if (opts.czt) {
	    fprintf(stderr, "# CZT mode\n");
	    fprintf(stderr, "# opts.n input windowed data size = %d\n", opts.n);
	    fprintf(stderr, "# opts.in input padded block size = %d\n", opts.in);
	    fprintf(stderr, "# opts.out output block size = %d\n", opts.out);
	    fprintf(stderr, "# opts.m output data size = %d\n", opts.m);
	    fprintf(stderr, "# output points to user = %d\n", opts.m);
	    fprintf(stderr, "# opts.singlesided = %d\n", opts.singlesided);
	    fprintf(stderr, "# opts.fflag=%d, opts.fstart=%g, opts.fstop=%g\n",
	    	opts.fflag, opts.fstart, opts.fstop);
	} else {
	    fprintf(stderr, "# FFT mode\n");
	    fprintf(stderr, "# opts.n input windowed data size = %d\n", opts.n);
	    fprintf(stderr, "# opts.in input padded block size = %d\n", opts.in);
	    fprintf(stderr, "# opts.out output points to user = %d\n", opts.out);
	    fprintf(stderr, "# opts.m output data size = %d\n", opts.m);
        }
        fprintf(stderr, "# window type = %s\n", name);
    }
}

double mag(COMPLEX z) {
   return(sqrt(z.re*z.re + z.im*z.im));
}

void putpoint(int i, double fcorrection, double foffset, double p) {
   double corr;
   if (opts.spect) {
      corr=sqrt(opts.wbw*opts.samplerate/(double) opts.n);
   } else {
      corr=1.0;
   }
   if (opts.dbmode) {
      printf("%12g %12g\n", 
      	foffset+((double)i)*fcorrection*opts.n/((double)opts.m-1.0), 
	20*log10(p*(opts.m/opts.n)/corr)-20*log10(opts.dbref));
   } else {
      printf("%12g %12g\n", 
      	foffset+((double)i)*fcorrection*opts.n/((double)opts.m-1.0), 
	p*(opts.m/opts.n)/corr);
   }
}

double max(double x, double y) {
    if (x>y) {
   	return x;
    } else {
       return y;
    }
}

int main(int argc, char **argv)
{
    extern OPTS opts;
    int i;
    COMPLEX *in;
    COMPLEX *spec;
    COMPLEX *out;
    COMPLEX *inbuf;
    COMPLEX *cztemp;
    double *win;
    double time, val;
    int nc;
    double fcorrection;
    double davg;
    double dc;		// returns dc value of each readin block
    int done = 0;
    int opt;
    int numwins = 0;
    double x;
    double gain;
    double m;
    int complex=0;	// set positive if readin() sees complex (three column) input

    init_opts();

    while ((opt = getopt(argc, argv, "01ab:d:f:hk:m:n:p:stw:v?")) != -1) {
	switch (opt) {
	case '0':
	    opts.nodc++;	// null dc component
	    break;
	case '1':
	    opts.singlesided++;	// plot singlesided spectrum
	    break;
	case 'a':		// animate: plot each avg as it is generated
	    opts.animate++;
	    break;
        case 'b':               // resolution bandwidth
	    opts.rbw = atof(optarg);
	    break;
	case 'f':		// frequency span
	    // FIXME: parse frequency span
	    opts.fflag = parse_fspec(optarg, &opts.fstart, &opts.fstop);
	    if (opts.fflag == -1) {
	       fprintf(stderr, 
	       	"error parsing -f frequency spec: %s\n", optarg);
	       exit(1);
	    }
	    opts.czt++;
	    break;
	case 'h':		// peak hold
	    opts.peakhold++;
	    break;
	case 'k':		// kaiser sidelobe depth
	    opts.kaiser = atof(optarg);
	    if (opts.kaiser < 15.0) {
	       fprintf(stderr, 
	       	"kaiser sidelobe <%s> too low, set to 15db\n", optarg);
		opts.kaiser = 15.0;
	    }
	    if (opts.kaiser > 180.0) {
	       fprintf(stderr, 
	       	"kaiser sidelobe <%s> too high, set to 180db\n", optarg);
		opts.kaiser = 180.0;
	    }
	    opts.wtype = wintype("kaiser");	// specify window type
	    win_set_sidelobe(opts.kaiser);	// and set sidelobe level
	    break;
	case 'm':		// chirp z output length
	    opts.czt++;
	    opts.m = atoi(optarg);
	    break;
	case 'n':		// input data block length
	    opts.n = atoi(optarg);
	    break;
	case 'p':		// pad factor
	    opts.p = atof(optarg);
	    break;
	case 's':		// spectral density mode
	    opts.spect++;
	    opts.pdplot++;	// turn on labels too
	    break;
	case 't':		// pdplot headers
	    opts.pdplot++;
	    break;
	case 'w':
	    opts.wtype = wintype(optarg);	// specify window type
	    win_set_sidelobe(opts.kaiser);	// and set sidelobe level
	    break;
	case 'd':		// db reference
	    opts.dbmode++;
	    opts.dbref = dbparse(optarg);	// specify db reference 
	    break;
	case 'v':		// verbose
	    printf("setting verbose\n");
	    opts.verbose++;
	    break;
	case '?':
	default:
	    fprintf(stderr, "%s: performs Welch's FFT averaging algorithm\n", argv[0]);
	    fprintf(stderr, "on standard input, writing spectrum plot to standard output.\n");
	    fprintf(stderr, "If -m<output block length> or -f<frequency range> options are\n");
	    fprintf(stderr, "specified, will optionally perform Chirp Z transform to allow\n");
	    fprintf(stderr, "efficient zooming into a desired frequency range.\n\n");
	    fprintf(stderr, "usage: %s [opts]  < stdin > stdout\n", argv[0]);
	    fprintf(stderr, "\t-0 (suppress DC component on each block)\n");
	    fprintf(stderr, "\t-1 (plot single-sided spectrum)\n");
	    fprintf(stderr, "\t-a (animate: plot each avg as generated)\n");
	    fprintf(stderr, "\t-b <resolution bw>\n");
	    fprintf(stderr, "\t-d<dbreference>\n");
	    fprintf(stderr, "\t\tcan be -dbm, -dbv, -dbc, -dbu, -db<n.n>\n");
	    fprintf(stderr, "\t-f [<fmin>,<fmax>|<fcenter>@<span>]\n");
	    fprintf(stderr, "\t-h (enable peak hold mode)\n");
	    fprintf(stderr, "\t-k <kaiser sidelobe width> default=120db\n");
	    fprintf(stderr, "\t-m <CZT output block length>\n");
	    fprintf(stderr, "\t-n <input data block length>\n");
	    fprintf(stderr, "\t-p <fft pad factor>\n");
	    fprintf(stderr, "\t-s (enable spectral density mode)\n");
	    fprintf(stderr, "\t-t (enable pdplot headers)\n");
	    fprintf(stderr, "\t-v (verbose mode)\n");
	    fprintf(stderr, "\t-w <window_type>\n");
	    fprintf(stderr, "\t\twindow can be (b)lackman, (han)ning, (ham)ming,\n");
	    fprintf(stderr, "\t\t(n)uttal, (k)aiser or (r)rectangular\n");
	    fprintf(stderr, "\t\tdefault is (b)lackman\n");
	    fprintf(stderr, "\t-? (print usage)\n");
	    exit(1);
	}
    }

    massage_opts();

    // Set up an array to hold the data, and assign the data.
    in = new_complex(opts.in);		// buffer for fft1d
    inbuf = new_complex(opts.in);	// scratch for readin();
    
    // preparing to change windowing technique
    // first create an array to hold window coefficients:
    win = new_double(opts.n);		// buffer for holding window
    //
    // then fill array with coefficients 
    //
    window_get(win,opts.n,opts.wtype);

    //for (i=0; i<opts.n; i++) {
    //   printf("# %g %g\n", (double)i, win[i]);
    //}

    // 
    // then get data as many times as we want in array in[]
    // window_do(in,win,opts.nodc,opts.singlesided);
    // (multiplies in[] by coefficients in win[])
    //


    if (opts.czt) {	// prepare czt auxiliary array
    	cztemp = new_complex(opts.out);
    	spec = new_complex(opts.m);	// hold avg power spect
    } else {
    	spec = new_complex(opts.in);	// hold avg power spect
    }

    // FIXME: readin() fails when number of input points is less than 
    // default 1024. should reset opts.n and try again.

    numwins=0;
    int nn=opts.n;
    // printf("#entering while %d\n", nn);
    while (1) {
        done=readin(&nn,inbuf,&davg,numwins,&complex);
	if (done && numwins != 0) break;
        
	// printf("#readin returned %d\n", nn);

	// if we run out of points on the first window, then adjust 
	// opts.n and continue

	if (numwins==0) {
	   if (nn != opts.n) {
	     fprintf(stderr,"ran out of points before reaching %d, continuing with %d\n", opts.n, nn);
	     opts.n=nn;
	   }
	}

	// don't overwrite original data for overlapping
	// windows:
        for (i = 0; i<opts.n; i++) {
	   in[i].re = inbuf[i].re; 
	   in[i].im = inbuf[i].im; 
	}
	numwins++;

	// dumparray(in, opts.n, "this is the original readin array");

        opts.samplerate = 1.0/davg;
	// fprintf(stderr,"samplerate is %g, davg=%g\n", opts.samplerate, davg);

	// old way:
	// creates window, corrects gain, every block...
	// window(in,opts.n,opts.wtype,opts.nodc,opts.singlesided);

	// new way:
	// creates window once in win[] and just reuses it...
        window_do(in,win,opts.n,opts.nodc,opts.singlesided);

	// dumparray(in, opts.n, "this is the windowed array");

	for (i = opts.n; i< opts.in; i++) {	// zero pad
	   in[i].re=0.0; in[i].im=0.0;
	}

	// dumparray(in, opts.in, "this is the zero padded array");

	if (!opts.czt) {		// normal fft
	    out=fft_1d(in,opts.in);
	    fcorrection = opts.samplerate/(double) opts.n;
	    // dumparray(in, opts.in, "this is the fftd array");
	} else {			// czt mode
	    m = (double) opts.m;
            switch (opts.fflag) {
		case 0:
		    // ((msamples-2)/msamples)*(-Fs/2) to Fs/2
		    // ((msamples-1)/msamples)*(-Fs/2)

		    opts.fstart = ((m-2.0)/(m))*(-opts.samplerate/2.0);
		    opts.fstop = opts.samplerate/(2.0);
		    break;
	    	case 1:	
		    opts.fstop = opts.samplerate*((m-1.0)/(2.0*m));
		    break;
		case 2:
		    opts.fstart = -opts.samplerate*((m-1.0)/(2.0*m));
		    break;
		case 3:
		    // FIXME: perhaps check for out of bounds frequencies here
		default:
		    break;
	    }
	    // fcorrection = opts.samplerate/m;
	    fcorrection = opts.p*(opts.fstop-opts.fstart)/m;

	    // fprintf(stderr,"fs=%g fe=%g fd=%g\n", opts.fstart, opts.fstop, opts.samplerate);
	    out=czt(in, 		//(COMPLEX*) input sequence
	            opts.in,    	//(int) size of input sequence
		    opts.m,     	//(int) size of output sequence
		    cztemp,     	//(COMPLEX*) array used to return results
		    opts.out,   	//(int) size of scratch array
		    opts.fstart,	//(double) start frequency
		    opts.fstop, 	//(double) stop frequency
		    opts.samplerate); 	//(double) sampling frequency
	}
	
	for (i = 0; i < opts.m; i++) {	
	    x = (mag(out[i])/(double)opts.m);
	    if (numwins == 1) {
		spec[i].re = x*x; 	// welch averages square of magnitude
		spec[i].im = 0.0;
	    } else {
	        if (opts.peakhold) { 		// save peaks
	       	    spec[i].re = max(spec[i].re, x*x); 
		    spec[i].im = 0.0;
	        } else {			// average sqr mags
		    spec[i].re += (x*x-spec[i].re)/((double)numwins);
		    spec[i].im = 0.0;
		}
	    }
	    if (opts.animate) {		// for now plots only half the spectrum
		putpoint(i, fcorrection, 0.0, sqrt(spec[i].re));
	    	// printf("  %12f  %12f\n", ((double) i)*fcorrection, spec[i].re);
	    }
	}
	// dumparray(spec, opts.m, "this is spec array");
    }

    if (opts.pdplot) {
	if (opts.dbmode) {
	    if (opts.spect) {
		printf("yscale 1 Magnitude (dBV/Hz)\n");
	    } else {
		printf("yscale 1 Magnitude dBv\n");
	    }
	} else {
	    if (opts.spect) {
		printf("yscale 1 Magnitude (Vrms/sqrt(Hz))\n");
	    } else {
		printf("yscale 1 Magnitude\n");
	    }
	}
        printf("xscale 1 Frequency   RBW=%g HZ  (%g dB)\n",
	        opts.wbw*opts.samplerate/(double) opts.n,
		10*log10(opts.wbw*fcorrection));
        // printf("title fs=%g rbw=%gHz, n=%d frames=%d\n",
        //	opts.samplerate, opts.wbw*fcorrection,opts.n,numwins);
    }

    if (!opts.czt && opts.singlesided && !complex) {
	double posmag, negmag, corr;
	for (i = 0; i < ((opts.m/2)+1); i++) {

	    // FOO
	    corr = 1.0; //sqrt(2.0);
	    corr = 1.0/sqrt(2.0);

	    if (i == 0) {		
		posmag = sqrt(spec[i].re);
		negmag = posmag;
	    } else {
		posmag = sqrt(spec[i].re);
		negmag = sqrt(spec[opts.m - i].re);
	    }

	    // if we are in nodc mode, the dc component is
	    // nulled out to about -350dB.  If we plot this
	    // point, it distorts the yaxis in decibal mode
	    // and makes it difficult to read the plot.  
	    // So if we are in dbmode, and nodc mode, we
	    // suppress the zero Hz component.

	    if (!(i==0 && opts.nodc && opts.dbmode)) {
	        putpoint(i, fcorrection, 0.0,
	    	    corr*(negmag + posmag));
	    }
	}
    } else {	// doublesided or czt
	if (!opts.czt) { // double sided
	    for (i = opts.m/2; i > 0; i--) {
		putpoint(-i, fcorrection, 0.0, sqrt(spec[opts.m - i].re));
	    }
	    // printf("pen\n");
	    for (i = 0; i < opts.m/2; i++) {
		putpoint(i, fcorrection, 0.0, sqrt(spec[i].re));
	    }
	} else { // czt
	    for (i = 0; i < opts.m; i++) {
		putpoint(i, fcorrection, opts.fstart, sqrt(spec[i].re));
	    }
	}
    }

    return(0);
}

//------------------------------------------------------------
// xfgets() is a replacement for fgets() with pushback.
// if xfgets_pushback(mode) is called with mode=1, then 
// lines are read from fp and saved in a linked list queue.
// if xfgets_pushback(0), the lines are read from the
// queue (in original order) until the queue is empty, and 
// then lines are read from fp again.

typedef struct strlist {
    char *s;
    struct strlist *next;
    struct strlist *prev;
} STRLIST;

static int pushback = 0;
static STRLIST *head = NULL;
static STRLIST *tail = NULL;

xfgets_pushback(int mode) {
   pushback=mode;
}

xfgets_dump() {
   STRLIST *p;
   for (p=head; p!=NULL; p=p->next) {
      fprintf(stderr,"%s\n", p->s);
   }
}

xfgets_save(char *s) {
   char *savebuf;
   STRLIST *strlist;
   STRLIST *tmp;

   // copy the string
   savebuf = malloc(strlen(s)+1);
   strcpy(savebuf, s);

   // create a STRLIST object
   strlist = malloc(sizeof(STRLIST));
   strlist -> next =  NULL;
   strlist -> prev =  NULL;
   strlist -> s = savebuf;

   if (head == NULL) {
      head = strlist;
      tail = strlist;
   } else {
      tail->next = strlist;
      strlist->prev = tail;
      tail = strlist;
   }
}

char *xfgets(char *buf, int size, FILE *fp) {
   char *s;
   STRLIST *tmp;

   if (pushback) {
       s=fgets(buf, MAXBUF, stdin);
       xfgets_save(s);
   } else {
       if (head != NULL) {
	   strcpy(buf,head->s);
           s = buf;
	   tmp = head;
	   head = head->next;
	   free(tmp);
       } else {
	   s=fgets(buf, MAXBUF, stdin);
       }
       return(s);
   }

   return(s);
}

//------------------------------------------------------------

// read in <n> doubles to array <in> return nonzero on error if eof
// computes running average of sample rate and returns
// average sample interval in <*davg>.  Gives error on
// stderr if any time spacing varies from the running average
// by more than TIMETOL.
// *complex is a pointer to a flag that reports it readin 
// ever seen a complex number on the input

#define TIMETOL 1e-1

int readin(int *nn, COMPLEX * in, double *davg, int numwins, int *complex) {
    int i;
    double t, told, delta,  re, im;
    int error=0;
    int done=0;
    static int line=0;
    char buf[MAXBUF];
    int retval;
    int offset;
    int n=*nn;

    offset=0;
    // printf("foo: in readin\n");
    if (numwins != 0) {	// not on first one, so shift and overlap
        for (i = 0; i <= n/2; i++) {
	   in[i].re=in[i+n/2].re;
	   in[i].im=in[i+n/2].im;
        } 
        offset=(n%2)+n/2;
    } else {
        offset=0;
    }


    for (i = 0; i < n-offset; i++) {
	told=t;

	// read complete lines and parses them
	// so that we accept 2 and three doubles per line
	// for complex data without slipping sync

	if (xfgets(buf, MAXBUF, stdin) == NULL) {
	    done++;
	    break;
	}
	line++;

	retval = sscanf(buf,"%lf %lf %lf", &t, &re, &im);

        if (retval==2) {
	    in[i+offset].re = re;
	    in[i+offset].im = 0.0;
	} else if (retval==3) {
	    in[i+offset].re = re;
	    in[i+offset].im = im;
	    *complex=1;
	} else {
	    fprintf(stderr, "line %d: bad fmt\n", line);
	    done++;
	    break;
	}

	delta=t-told;
	if (i==2) {
	   *davg = delta;
	} else {
           *davg += (delta-*davg)/((double)i);
	}

	// printf("## %g %g\n",delta, *davg);

	if (error == 0 && 
		n > 2 && 
		fabs((delta-*davg)/(delta+*davg)) > TIMETOL) {
            error++;		// suppress second error
     	    fprintf(stderr,
	    	"# WARNING: bad timestep at line %d, %g %g %g %g\n",
		line, t, told, delta, *davg);
	}  else { 
	    error=0;
        }
    }
    // printf("foo: exiting readin with %d\n",i);
    if (numwins==0) *nn = i;
    return(done);
}

/*
        -d ref By default, the spectrum is plotted in rms magnitude,
        with the pretense that the signal has the dimension of Volts.
        The -d option directs welch to output the power spectrum in
        decibels.  The argument ref can be the strings "bv", "bm", "bu",
        or "bc", or a number may be given, specifying output in dBV
        (rms), dBm (into 50 ohms), dBu (standard in audio electronics --
        Vref=sqrt(0.6) Vrms, the voltage that would dissi‐ pate 1mW into
        600 ohms), dBc (dB relative to the output signal's peak
        magnitude), or dB (rms) relative to the given number.  For
        example, "-dbm" plots the power in dBm and "-d1e-3" plots the
        magnitude in decibels relative to 1mV.
*/


double dbparse(char *dbrefstring) {
    double dbref;
    if (strncasecmp(dbrefstring,"bv",3) == 0) {
        dbref = 1.0;
    } else if (strncasecmp(dbrefstring,"bm",3) == 0) {
        dbref = 50e-3;
    } else if (strncasecmp(dbrefstring,"bu",3) == 0) {
        dbref = 600e-3;
    } else {
        fprintf(stderr,"# bad db reference \"d%s\", using dbv\n", optarg);
        dbref = 1.0;
    }
    // printf("called with %s, return %d\n",wname, type);
    return(dbref);
}

int parse_fspec(char *line, double *fstart, double *fstop) {
   int i;
   int nval;
   char buf[128];
   char *str1 = NULL;
   char *str2 = NULL;
   double tmp;
   int debug=0;
   int retcode = 0;

   if (debug) fprintf(stderr,"testing %s\n", line);
   strcpy(buf,line);

   str1=buf; str2=NULL;
   for (i=0; buf[i]!='\0' && i<(strlen(buf)-1); i++) {
      if ((buf[i]==',') && (buf[i+1]!='\0')) {
	 buf[i]='\0';	// null terminate first line
	 str2=buf+i+1;	// set str2 to second line
	 break;
      }
   }

   if (debug) fprintf(stderr,"str1=<%s>, str2=<%s>\n", str1, str2);

   if (str1 != NULL && strlen(str1) > 0) {
       nval=sscanf(str1, "%lf", &tmp);
       if (nval==1) {
          *fstart = tmp; 
	  retcode+=1;
       } else {
	  return(-1);	// unparseable 
       }
   }
   if (str2 != NULL && strlen(str2) > 0) {
       nval=sscanf(str2, "%lf", &tmp);
       if (nval==1) {
          *fstop = tmp; 
	  retcode+=2;
       } else {
	  return(-1);	// unparseable 
       }
   }
   if (debug) fprintf(stderr,"code= %d, fstart=<%g>, fstop=<%g>\n", retcode, *fstart, *fstop);
   return (retcode);
}

void dumparray(COMPLEX *in, int npts, char *str) {
   int i;
   fprintf(stderr,"# %s\n", str);
   fprintf(stderr,"# npts=%d\n",npts);
   for (i=0; i<npts; i++) {
      fprintf(stderr,"%d %g %g\n", i, in[i].re, in[i].im);
   }
}
