
typedef struct cplex {
   double re;
   double im;
} COMPLEX;

extern COMPLEX *new_complex(int size);
extern COMPLEX *fft_1d( COMPLEX *array, int n);
extern COMPLEX *ifft_1d(COMPLEX *array, int n);

extern COMPLEX *czt( 
	COMPLEX *array, 
	int N, 
	int M, 
	COMPLEX *g,
	int L,
	double fstart, 
	double fstop, 
	double fsam
    );

extern int snap2(int n);  // snap to next enclosing power of two
extern int czt_blocksize(int N, int M);
