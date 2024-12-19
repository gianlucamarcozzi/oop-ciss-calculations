#define	PI	M_PI
#define DataT double
#define MAX_NUMBER_POINTS 32
#define BE 1.399624167     /* beta/h; (vorher 1.39968)  SZ 9/8/96 */
#define GBETA (2.0023193*BE) 
#define FACTOR 0.398942   /*(1/sqrt(2*PI))*/
#define SQR(x)	((x)*(x))
#define asc ".asc\0"
#define deg (180.0/PI)
#define rad (PI/180.0)
/*#define register  short*/  /* SZ 26/8/95 */

typedef DataT
   vector[3],
   vector4[4],
   vector10[10];

typedef vector
   matrix[3];

typedef char
   strg[16];

strg		filename[9], orient_name[3];

double   	*spec, **Transpec, **tens, *orient_spec;

double		*tensor;

matrix      	aA1, aP, gP, gP2, gA1, gA12, gaP2, ga2, Dm,
                rPQ, rtPQ, rCQ, rtCQ, rBC, rtBC, Ra,
		Rz180, Rz120, Rz240;

vector		rpA1, euler_PQ, euler_zP, euler_zQ, zP, zQ;
                

vector4		intens, trans;

vector10	hyperfine_prob_P, hyperfine_prob_A1;

DataT		fmw, h1, h2, dh, phi1, theta1, J, d, Jpd, Jmd, geffP,
		geffA1, aeffP, aeffA1, hbP, hbA1, hbPA1, miA1, miP,
		lorentz_hb, max_freq;

int 		imiA1, imiP, nt, np, step_size, number_of_protons,
		n_protons_p, NumOfPoints, delta_w_steps, num_freq_points,
		alpha_step, astep, num_spec;

char    	reference_syst, print_orient, print_Powder,
		print_ESEEM_FT,  outfile[20];
				
    
