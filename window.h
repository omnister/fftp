
typedef enum {BLACK, HANN, HAMM, RECT, TRIANGULAR, NUTT, BNUT, KAISER, FLAT, NONE} WINTYPE;

extern int wintype(char *wname);				// return WINTYPE from a string
extern int win_set_sidelobe(double sidelobe);			// set kaiser sidelobe level
extern int wintype_gain(WINTYPE wintype, double *gain);		// get gain for current window
extern int wintype_rbw(WINTYPE wintype, double *rbw);		// get rbw for current window
extern int wintypetoname(WINTYPE wintype, char **wname);	// get canononic window name
extern int window(COMPLEX *in, int N, WINTYPE type, int nulldc, int singlesided);  	// window the *in array

