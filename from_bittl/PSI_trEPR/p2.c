/***************************************************************
*            p2.c                                             *
*          cc p2.c -o p2 -lm                                 *
****************************************************************
*  Versions:
*       - ESEEM_FTs and Tr EPR Spectra! (21.08.98 SZ)
*       - ESEEM with or without Hfc's (27.08.98 SZ)        
*       - Whole program derived vom p2k.c (version 31.08.98)
*
****************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "p2.h"
#include "p2_utils.c"


int ReadInput()
{
   short   i, k;
   FILE       *fp, *fopen();

   if ((fp = fopen("my_p2.txt","r")) == NULL) {
      printf("file p2.txt cannot be opened");
      return(-1);
   }
   fscanf(fp,"%lf %*s %lf %lf %d %*s %d %d %*s", &fmw, &h1, &h2, 
                                           &NumOfPoints, &nt, &np);
   dh = (h2-h1)/(NumOfPoints-1);
   for (i=0;i<=2;i++)
     for (k=0;k<=2;k++)
       fscanf(fp, "%lf ", &gP[i][k]);
   NewLine(fp);
   fscanf(fp,"%lf", &hbP);
   NewLine(fp);
   for(i=0;i<=2;i++)
     for (k=0;k<=2;k++)
       fscanf(fp, "%lf", &gA1[i][k]);
   NewLine(fp);
   fscanf(fp, "%lf", &hbA1);
   NewLine(fp);
   reference_syst=getc(fp);
   NewLine(fp);
   switch (toupper(reference_syst)) {
	case 'Q':
	  fscanf(fp, "%lf %lf %*s",&phi1,&theta1);
	  fscanf(fp, "%lf %lf %lf %*s",&euler_PQ[0],&euler_PQ[1],&euler_PQ[2]);
	  break;
/*
	case 'Z':
	  fscanf(fp, "%lf %lf %lf %*s",&euler_zP[0],&euler_zP[1],&euler_zP[2]);
	  fscanf(fp, "%lf %lf %lf %*s",&euler_zQ[0],&euler_zQ[1],&euler_zQ[2]);
	  break;
	case 'P':
	  fscanf(fp, "%lf %lf %*s",&phi1,&theta1);
	  fscanf(fp, "%lf %lf %lf %*s",&euler_PQ[0],&euler_PQ[1],&euler_PQ[2]);
	  break;
*/
	default:
	  printf("Up to now, only Q as reference system possible (SORRY!)  ");
	  fclose(fp);
	  return -1;
   }
   fscanf(fp, "%lf %lf %*s", &J, &d);
   NewLine(fp);
   print_orient = getc(fp);
   NewLine(fp);
   print_Powder = getc(fp);
   NewLine(fp);
   print_ESEEM_FT = getc(fp);
   NewLine(fp);
   fscanf(fp, "%20s %*s", outfile);
   fscanf(fp, "%lf %d %*s", &max_freq, &num_freq_points);
   fscanf(fp, "%d %*s", &delta_w_steps);
   fscanf(fp, "%lf %*s", &lorentz_hb);
   fscanf(fp,"%d %*s",&number_of_protons);
   for(i=0;i<=2;i++)
     for (k=0;k<=2;k++)
       fscanf(fp, "%lf", &aA1[i][k]);
   fscanf(fp,"%d %*s",&n_protons_p);
   for(i=0;i<=2;i++)
     for (k=0;k<=2;k++)
       fscanf(fp, "%lf", &aP[i][k]);
   fclose(fp);
   return(0);
}


void WriteInput()
{
    short 	i, k;

    printf("CCRP simulation for P+Q- powder spectra\n");
    printf("Microwave frequency: %.2f MHz \n",fmw);
    printf("Magnetic field: from %.2f to %.2f Gauss in steps of %.4f Gauss\n"
                                                    ,h1,h2,dh);
    printf("\nNumber of Fieldpoints = %d \n",NumOfPoints);
    printf("Theta steps: %d      Phi steps: %d\n",nt,np);
    printf("\n ANGLES:\n");
    printf("zPQ_gQ:     Phi1: %f      Theta1: %f\n ",phi1, theta1);
    printf("Euler(gP_gQ) %f, %f, %f\n ",euler_PQ[0],euler_PQ[1], euler_PQ[2]);
    PrintMatrix(gP,"g-Tensor for P+");
    PrintMatrix(gA1,"g-Tensor for A1-");
    printf("Linewidth P: %.2f  Q: %.2f\n",hbP,hbA1);
    printf("Dipol-Dipol-Coupling %.3f Exchange Coupling %.3f\n", d,J);
    if (number_of_protons != 0) {
	printf("Acceptor Q has %d aequivalent protons\n",number_of_protons);
	PrintMatrix(aA1,"Hyperfine tensor");
    }
    if (n_protons_p != 0) {
	printf("Donor P has %d aequivalent protons\n",n_protons_p);
	PrintMatrix(aP,"Hyperfine tensor");
    }
    printf("Output filename %-20s\n",outfile);
}

void  PrintOutRelOrntPQ(o_sys)
     char o_sys;
{
  printf("\n RELATIVE ORIENTATION PQ \n");
  PrintMatrix(rPQ,"Transformation matrix P->Q");
  PrintMatrix(rtPQ,"Transformation matrix Q->P");
  printf("in g%c system:\n",reference_syst);
  printf("                alpha = %.2f\n",euler_PQ[0]*deg);
  printf("                beta  = %.2f\n",euler_PQ[1]*deg);
  printf("                gamma = %.2f\n",euler_PQ[2]*deg);
  printf("in g%c system:\n",o_sys);
  printf("                alpha = %.2f \n",180-euler_PQ[2]*deg);
  printf("                beta  = %.2f \n",euler_PQ[1]*deg);
  printf("                gamma = %.2f \n",180-euler_PQ[0]*deg);
}

void PrintOutAngles_zPQ(sys, theta, phi, vec)
     char     sys;
     DataT    theta, phi;
     vector   vec;
{
  char  axis[9];

  switch (toupper(sys)){
  case 'Q':
    strncpy(axis, "zPQ", sizeof(axis));
    printf("\n ORIENTATION of Dipolar Vector relative to g%c\n",sys);
    break;
  case 'P':
    strncpy(axis, "zPQ", sizeof(axis));
    printf("\n ORIENTATION of Dipolar Vector relative to g%c\n",sys);
    break;
/*
  case 'C':
    strncpy(axis, "c-axis", sizeof(axis));
    sys = 'Q';
    printf("\n ORIENTATION of crystall c-axis relative to g%c\n",sys);
    break;
  case 'D':
    strncpy(axis, "c-axis", sizeof(axis));
    sys = 'P';
    printf("\n ORIENTATION of crystall c-axis relative to g%c\n",sys);
    break;
*/
  case 'B':
    strncpy(axis, "B0-field", sizeof(axis));
    sys = 'Q';
    printf("\n ORIENTATION of B0-field relative to g%c\n",sys);
    break;
  default:
    return;
  }
  printf("  theta = %.2f (%.2f)  phi = %.2f (%.2f) \n",
               theta*deg,180-theta*deg,phi*deg,180-phi*deg);
  printf("  Angle: %s g%cx = %.2f (%.2f)\n", axis, sys,
                vec[0]*deg,180-vec[0]*deg);
  printf("  Angle: %s g%cy = %.2f (%.2f)\n", axis, sys,
                vec[1]*deg,180-vec[1]*deg);
  printf("  Angle: %s g%cz = %.2f (%.2f)\n", axis, sys,
                vec[2]*deg,180-vec[2]*deg);
}


/* calcu. rotation matrix from Euler anlges 
   R = Rz(gamma) Ry(beta) Rz(alpha)         */
void GetR(r,euler)      
     matrix     r;
     vector  euler;
{
  matrix   a, b, c;

  RotationMatrix(a,euler[0]*deg,'z');
  RotationMatrix(b,euler[1]*deg,'y');
  RotationMatrix(c,euler[2]*deg,'z');
  MulMat(r,c,b);
  MulMat(r,r,a);
}

/* Dipolar axis in Q axis system 
  rpA1 = Rz(phi) Ry(theta) (0,0,1)   */
void  GetrpA1()
{
    rpA1[0] = -sin(theta1) * cos(phi1);
    rpA1[1] =  sin(theta1) * sin(phi1);
    rpA1[2] =  cos(theta1);
}

/* Angle between dipolar axis and Q axes system */
void GetzQ()
{
    zQ[0] =  acos(rpA1[0]);
    zQ[1] =  acos(rpA1[1]);
    zQ[2] =  theta1;
}

/* Angle between dipolar axis and P axes system */
void GetzP()
{
    matrix	r;
    short	i,k;

    GetR(r,euler_PQ);
    Transpose(r,r);
    for(i=0;i<3;i++) {
       for(zP[i]=0.0, k=0; k<3; k++)
	   zP[i] += r[i][k]*rpA1[k];
       zP[i] = acos(zP[i]);
    }
}

void GetAngles()
{
  short  i;
  char      other_sys;
  float     phiP, phiC;

  if (reference_syst == 'P')
    other_sys = 'Q';
  else
    other_sys = 'P';
  GetrpA1();                     /* ----> zPQ in gQ axis system */
  GetzQ();
  GetzP(); 
  PrintOutAngles_zPQ('Q',theta1,phi1,zQ);
  phiP = atan(-cos(zP[1])/cos(zP[0]));
  PrintOutAngles_zPQ('P',zP[2],phiP,zP);
  PrintOutRelOrntPQ(other_sys);
}


/* Transf. gP in gQ axis system */
void RotateGP()   
{
  PrintMatrix(gP," gP in gP");
  GetR(rPQ, euler_PQ);
  Transpose(rtPQ, rPQ);
  MulMat(gP,gP,rtPQ);
  MulMat(gP,rPQ,gP);
  PrintMatrix(gP," gP in gQ");
}


void GetHyperfineProbs(hyperfine_prob,np)
vector10  hyperfine_prob;
int 	  np;
{
   vector10	b;
   double	hnorm;
   short	i, j;

   for(i=0; i<=np+1; i++)
      hyperfine_prob[i] = 0.0;
      hyperfine_prob[1] = 1.0;
   if (np != 0 )  {
      for (b[0] = 0,i=1; i<=np+1; i++) {
	  for (j=1; j<=i; j++)
	      b[j] = hyperfine_prob[j] + hyperfine_prob[j-1];
	  for (j=1; j<=i; j++)
	      hyperfine_prob[j] = b[j];
      }
      for (hnorm=0,i=0;i<10;i++)
	  hnorm += hyperfine_prob[i];
      for (i=0;i<10;i++)
	  hyperfine_prob[i] /= hnorm;
   }
}

void RotDipolarMatrix()
{
   short i, j;
   matrix Rd, Rdt, Rdz, Rdy; /*, Djunk1, Djunk2; */
  
   for (i=0; i<3; i++)
      for(j=0; j<3; j++)
         Dm[i][j] = 0.0;
   Dm[0][0] = -d/3.0;
   Dm[1][1] = -d/3.0;
   Dm[2][2] = 2*d/3.0;
   PrintMatrix(Dm, "Dipol-Tensor vor Drehung");
   /* Rotationsmatrix: 1. Rot. um y-Achse um theta, 2. Rot. um z um phi 
          D' = Rz(phi) Ry(theta) D Ryt(theta) Rzt(phi) */
   RotationMatrix(Rdy, theta1 * deg, 'y');
   RotationMatrix(Rdz, phi1 * deg, 'z');
   MulMat(Rd, Rdz, Rdy);
   Transpose(Rdt,Rd);
   MulMat(Dm, Rd, Dm);
   MulMat(Dm, Dm, Rdt);
/* vandenBrink et al. 1995 wieder gel�scht. SZ 20/08/98
   for(i=0; i<3; i++)
      for(j=0; j<3; j++)
	 Djunk1[i][j] = 0.0;
   MulMat(Djunk1, gA12, Dm);
   MulMat(Djunk1, Djunk1, gP2);
   for(i=0; i<3; i++)
      for(j=0; j<3; j++)
	 Djunk2[i][j] = 0.0;
   MulMat(Djunk2, gP2, Dm);
   MulMat(Djunk2, Djunk2, gA12);
   for (i=0; i<3; i++)
      for(j=0; j<3; j++)
	 Dm[i][j] = (Djunk1[i][j] + Djunk2[i][j])/2.0;
   PrintMatrix(Dm, "Dipol-Tensor 1/2 (gA2 D gP2 + gP2 D gA2");
   MulMat(Dm, gA12, Dm);
   MulMat(Dm, Dm, gP2);
   PrintMatrix(Dm, "Dipol-Tensor nach g-Mat. mult.");
*/
}   

/* Calc. actual dipolar coupling: deff = Bvec  D Bvec  */ 
double Get_deff(B0vector)
vector B0vector;
{
   short i, j;
   double deff;              /*, teiler; */

   /* teiler = geffA1*geffP*4.0092053; */
   for (deff=0.0, i=0; i<3; i++) 
       for (j=0; j<3; j++)
		  deff += B0vector[i] * B0vector[j] * Dm[i][j];
   return(deff);
   /* Bei Multiplikation mit G-tensoren: Normierung  
   return(deff/teiler); 
   */
}



void InitPar()
{
   double  ct, st, sp, cp;
   short  i;
   matrix gPt, gA1t;

   step_size = nt*np;
   J = GBETA*J;
   d = GBETA*d;
   hbP = 0.5*SQR(GBETA*hbP);
   hbA1 = 0.5*SQR(GBETA*hbA1);
   hbPA1 = (hbP + hbA1)/2.0;
   theta1 = theta1 * rad;
   phi1 = phi1 * rad;
   for (i=0; i<3; i++){
       euler_PQ[i] *= rad;
   }
   RotateGP();
   GetAngles();
   Transpose(gPt, gP);
   Transpose(gA1t, gA1);
   MulMat(gP2,gP,gPt);
   MulMat(gA12,gA1,gA1t);
   RotDipolarMatrix();
   PrintMatrix(gP2,"g-Tensor fuer P2+");
   PrintMatrix(gA12,"g-Tensor fuer A12");
   if (number_of_protons != 0) {
      MulMat(ga2,gA1,aA1);
      MulMat(ga2,ga2,ga2);
   }
   if (n_protons_p != 0) {
      MulMat(gaP2,gP,aP);
      MulMat(gaP2,gaP2,gaP2);
   }
   GetHyperfineProbs(hyperfine_prob_P,n_protons_p);
   PrintVec10(hyperfine_prob_P," Hyperfein-Wahrscheinlichkeiten fuer P");
   GetHyperfineProbs(hyperfine_prob_A1,number_of_protons);
   PrintVec10(hyperfine_prob_A1," Hyperfein-Wahrscheinlichkeiten fuer A1");
}

void Getgeff(B0vector)
vector B0vector;
{
   double deff;
   short  i, k;

   for (geffA1=geffP=0,i=0; i<3; i++) {
       for (k=0; k<3; k++)
           geffA1 += B0vector[i] * B0vector[k] * gA12[i][k];
       for (k=0; k<3; k++) 
	   geffP += B0vector[i] * B0vector[k] * gP2[i][k];
   }
   geffP = sqrt(geffP);
   geffA1 = sqrt(geffA1);
   if (number_of_protons != 0) {
      for (aeffA1=0,i=0; i<3; i++)
	  for (k=0; k<3; k++)
	     aeffA1 += B0vector[i] * B0vector[k] * ga2[i][k];
      aeffA1 = sqrt(aeffA1)/geffA1;
   }
   if (n_protons_p != 0) {
      for (aeffP=0,i=0; i<3; i++)
	  for (k=0; k<3; k++)
	     aeffP += B0vector[i] * B0vector[k] * gaP2[i][k];
      aeffP = sqrt(aeffP)/geffP;
   }
   deff = Get_deff(B0vector); 
   Jpd = J + deff/2.0;
   Jmd = J - deff;
}


void  GetTransitionFrequencies(h)
double 	h;
{
   double w0, q, inten, p_character, hb12, hb34, norm12, norm34, x;

   w0 = BE * h * (geffA1 + geffP) / 2.0
	     + miA1 * aeffA1/2.0 + miP * aeffP/2.0;
   q =  BE * h * (geffP-geffA1)/2.0
	    -miA1 * aeffA1/2.0 + miP * aeffP/2.0;

   x = sqrt(Jpd*Jpd + q*q)*q/fabs(q); 
   p_character = 0.5 * (1 + q/x);
   hb12 = p_character * hbP + (1. - p_character) * hbA1;
   hb34 = p_character * hbA1 + (1. - p_character) * hbP;
   trans[0] = -SQR(-fmw + w0 + Jmd + x)/hb12;
   trans[1] = -SQR(-fmw + w0 - Jmd + x)/hb12;
   trans[2] = -SQR(-fmw + w0 + Jmd - x)/hb34;
   trans[3] = -SQR(-fmw + w0 - Jmd - x)/hb34;
   inten =  q*q/(4 * (Jpd*Jpd + q*q))
		 * hyperfine_prob_A1[imiA1] * hyperfine_prob_P[imiP];
   intens[0] = -inten * FACTOR / hb12;
   intens[1] =  inten * FACTOR / hb12;
   intens[2] = -inten * FACTOR / hb34;
   intens[3] =  inten * FACTOR / hb34;
}


void SingleOrientationSpectrum(B0vector)
vector B0vector;
{
   double     h;
   short   field_pt, k;

   Getgeff(B0vector);
   for (miA1 = -number_of_protons/2.0,imiA1 = 1;
	     imiA1 <= number_of_protons + 1;  imiA1++) {
      for (miP = -n_protons_p/2.0, imiP = 1;
	       imiP <= n_protons_p + 1;  imiP++) {
	 for (h=h1,field_pt=0; field_pt < NumOfPoints; field_pt++) {
	    GetTransitionFrequencies(h);
	    if ((imiA1 == 1) && (imiP == 1))
	       orient_spec[field_pt] = 0.0;
	    for (k = 0; k < 4; k++) {
	       if (trans[k] > -15.0)  {
		  Transpec[k][field_pt] += intens[k]* exp(trans[k])/step_size;
		  orient_spec[field_pt] += intens[k]* exp(trans[k])/step_size;
	       }
	    }
	    if ((imiA1 == number_of_protons + 1) && (imiP == n_protons_p +1))
	       spec[field_pt] += orient_spec[field_pt];
	    h += dh;
	 }
	 miP += 1.0;
      }
      miA1 += 1.0;
   }
}



/* Direction of B0 in Q axes system */
void GetB0vector(B0vector, cos_theta_B0, phi_B0)
     vector B0vector;
     double cos_theta_B0, phi_B0;
{
  double sin_theta_B0;

  sin_theta_B0 = Compliment(cos_theta_B0);
  B0vector[0] = sin_theta_B0 * cos(phi_B0);
  B0vector[1] = sin_theta_B0 * sin(phi_B0);
  B0vector[2] = cos_theta_B0;
}


void MakeTensorWithHfc(B0vector)
vector B0vector;
{
    double  freq, faktor, A, B, breite, delta_w, delta_w_i, R_2, R;
    double exponent, y, inten;
    long int  freq_point, i;

    Getgeff(B0vector);
    A = 2*Jmd;
    B = 2*Jpd;
    delta_w = BE * 0.5*(h1+h2) * (geffP - geffA1);
    breite = sqrt(2*hbPA1);
    for (delta_w_i = delta_w - 3*breite ; delta_w_i <= delta_w + 3*breite;
		delta_w_i += 6*breite/delta_w_steps) {
        exponent = -SQR(delta_w - delta_w_i)/hbPA1;  
	inten = exp(exponent);
	R_2 = SQR(delta_w_i) + SQR(B);
        if (delta_w_i != delta_w)
	   R = sqrt(R_2) * (delta_w_i-delta_w)/fabs(delta_w_i-delta_w); 
        else
           R = sqrt(R_2);
	if (R_2 != 0.0)
	   faktor = SQR(delta_w_i) * SQR(B) / SQR(R_2);
        for (freq_point = 0, freq = - max_freq ; 
	                       freq_point < num_freq_points; freq_point++) {    
	   if ((y = SQR((A-freq)/lorentz_hb)) < 1e+03) {
	      tens[0][freq_point] += faktor/(1 + y);
	      tensor[freq_point] += faktor/(1 + y);
           }
	   if ((y = SQR((A-R-freq)/lorentz_hb)) < 1e+03) {
	      tens[1][freq_point] += -0.5*inten*faktor/(1 + y);
	      tensor[freq_point] += - 0.5*inten*faktor/(1 + y);
           }
	   if ((y = SQR((A+R-freq)/lorentz_hb)) < 1e+03) {
	      tens[2][freq_point] += -0.5*inten*faktor/(1 + y);
	      tensor[freq_point] += -0.5*inten*faktor/(1 + y);
	   }

	freq += 2*max_freq/(num_freq_points-1);
	}
    }
}

/* ESEEM FTs without consideration of Hfc's ot G-anisotropy --
   only dipolar and exchange interaction  Gamma = 2 (J - d) 
   with Lorentzian broadening */

void MakeTensor(B0vector)
vector B0vector;
{
    double  freq, Gamma, y; 
    long int  freq_point, i;

    Getgeff(B0vector);
    Gamma = 2*Jmd;

    for (freq_point = 0, freq = - max_freq ; freq_point < num_freq_points; 
                                                      freq_point++) {
           if ((y = SQR((Gamma-freq)/lorentz_hb)) < 1e+03) {
              tensor[freq_point] += 1/(1 + y);
           }
        freq += 2*max_freq/(num_freq_points-1);
    }
}

/* Add Gamma and -Gamma contribution (Point Symmetry of SFT) */
void MakeESEEMFourier()
{
   double *temp;
   int i, k;

   if ((temp = (double *)malloc(num_freq_points*sizeof(double))) == NULL) {
      printf("Couldn't allocate memory for TEMP\n");
   }

   for (i=0; i<num_freq_points; i++)
      temp[i] = tensor[num_freq_points - 1 - i];
   for (i=0; i<num_freq_points; i++) 
      tensor[i] -= temp[i];
}

void MakeESEEMFourierWithHfc()
{
   double *temp;
   int i, k;

   if ((temp = (double *)malloc(num_freq_points*sizeof(double))) == NULL) {
      printf("Couldn't allocate memory for TEMP\n");
   }

   for (i=0; i<num_freq_points; i++) 
      temp[i] = tensor[num_freq_points - 1 - i];
   for (i=0; i<num_freq_points; i++) { 
      tensor[i] -= temp[i];
      /*
      tensor[i] /= delta_w_steps;
      */
   }
   for (k=0; k<3; k++)  {
      for (i=0; i<num_freq_points; i++)
	 temp[i] = tens[k][num_freq_points - 1 - i];
      for (i=0; i<num_freq_points; i++) {
	 tens[k][i] -= temp[i];
	 /*
         tens[k][i] /= delta_w_steps;
	 */
      }
   }
}       

void Write_FT_Files()
{
    char spec_fn[20];

       sprintf(spec_fn, "%s_ten", outfile);
       printf("%s\n", spec_fn); 
       WriteSpec(-max_freq, max_freq, tensor, num_freq_points, spec_fn); 
       sprintf(spec_fn, "%s_A", outfile);
       WriteSpec(-max_freq, max_freq, tens[0], num_freq_points, spec_fn);
       sprintf(spec_fn, "%s_AmR", outfile);
       WriteSpec(-max_freq, max_freq, tens[1], num_freq_points, spec_fn);
       sprintf(spec_fn, "%s_ApR", outfile);
       WriteSpec(-max_freq, max_freq, tens[2], num_freq_points, spec_fn);
}    

void Write_TrEPR_Files()
{
    int i,k;
    char spec_fn[20];

       sprintf(spec_fn, "%s_spec", outfile);
       printf("%s\n", spec_fn);
       WriteSpec(h1, h2, spec, NumOfPoints, spec_fn);
       sprintf(spec_fn, "%s_Pem", outfile);
       WriteSpec(h1, h2, Transpec[0], NumOfPoints, spec_fn);
       sprintf(spec_fn, "%s_Pabs", outfile);
       WriteSpec(h1, h2, Transpec[1], NumOfPoints, spec_fn);
       sprintf(spec_fn, "%s_Qem", outfile);
       WriteSpec(h1, h2, Transpec[2], NumOfPoints, spec_fn);
       sprintf(spec_fn, "%s_Qabs", outfile);
       WriteSpec(h1, h2, Transpec[3], NumOfPoints, spec_fn);
       for (i = 0; i < NumOfPoints; i++) {
	   spec[i] = 0.;
	   for (k = 0; k < 4; k++) spec[i] += Transpec[k][i];
       }
       sprintf(spec_fn, "%s_sum", outfile);
       WriteSpec(h1, h2, spec, NumOfPoints, spec_fn);
}

void MakeSpec()
{
  double     cos_theta_B0, phi_B0, delta_cos_theta, delta_phi;
  vector     B0vector;
  short      i, i_cos_theta, i_phi;

  delta_cos_theta = 2.0/nt;
  delta_phi = PI/np;

  for (cos_theta_B0 = 1 - delta_cos_theta/2, i_cos_theta = 1;
                 i_cos_theta <= nt; i_cos_theta++) {
     printf("%d\t%f\n", nt - i_cos_theta, cos_theta_B0); fflush(stdout);	  
     for (phi_B0 = 0.0, i_phi=1; i_phi <= np; i_phi++) {
        GetB0vector(B0vector, cos_theta_B0, phi_B0);       
	if (tolower(print_Powder) == 'y')
           SingleOrientationSpectrum(B0vector);       
        if (tolower(print_ESEEM_FT) == 'y') {
           if (delta_w_steps == 0) {
              MakeTensor(B0vector);
           }   
           else {
              MakeTensorWithHfc(B0vector); 
	   }
        }
	phi_B0 += delta_phi;      
    }
    cos_theta_B0 -= delta_cos_theta;
  }  
  if (tolower(print_Powder) == 'y')
     Write_TrEPR_Files();
  if (tolower(print_ESEEM_FT) == 'y') {
     if (delta_w_steps == 0)
        MakeESEEMFourier();
     else 
        MakeESEEMFourierWithHfc();
     Write_FT_Files();
  }
}

void RingBell(n)
int n;
{
   short i;

   for (i=1; i<=n; i++) {
      printf("\7 \n");
   }
   printf("Fertig! \n");
}


void InitVariables(NumOfPoints, FreqPkt)
int NumOfPoints, FreqPkt;
{
   int i, k;


   if ((tensor = (double *)malloc(FreqPkt*sizeof(double))) == NULL) 
      printf("Couldn't allocate memory for TENSOR\n");

   if ((tens = (void *)malloc(3*sizeof(double *))) == NULL)
      printf("No Memory for TENS-Pointer\n");
      
   for (k=0; k<3; k++) {
      if ((tens[k] = (double *)malloc(FreqPkt*sizeof(double))) == NULL)
	 printf("Couldn't allocate memory for Tens\n");
      for (i=0; i<FreqPkt; i++)
	 tens[k][i] = 0.0;
   }

   for (i = 0; i < FreqPkt; i++)  
      tensor[i] = 0.0;

   if ((Transpec = (void *)malloc(4*sizeof(double *))) == NULL)
      printf("No Memory for Transpec-Pointer\n");


   for (k=0; k < 4; k++) {
      if ((Transpec[k] = (double *)malloc(NumOfPoints*sizeof(double))) == NULL) 
	 printf("Couldn't allocate memory for Transpec\n");
      for (i=0; i < NumOfPoints; i++) 
	 Transpec[k][i] = 0.0;
   }

   if ((spec = (double *)malloc(NumOfPoints*sizeof(double))) == NULL)
      printf("Couldn't allocate memory for SPEC\n");
   if ((orient_spec = (double *)malloc(NumOfPoints*sizeof(double))) == NULL)
      printf("Couldn't allocate memory for ORIENT_SPEC\n");
   for (i=0; i < NumOfPoints; i++) {
      spec[i] = 0.0;
      orient_spec[i] = 0.0;
   }
}

int main()  {

   int k; 

   if (ReadInput() == 0) {
      WriteInput();
      InitVariables(NumOfPoints, num_freq_points);
      InitPar();
      MakeSpec();  
      /*
      WriteGnuInput();
      */
      RingBell(1);
   }
}

