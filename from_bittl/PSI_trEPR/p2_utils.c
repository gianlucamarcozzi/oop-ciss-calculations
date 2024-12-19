void PrintMatrix(a,title)
matrix	a;
char	title[];
{
   int	i, j;

   printf("\n         %s \n",title);
   for (i=0;i<=2;i++) {
      printf("     ");
      for (j=0;j<=2;j++)
	 printf("%9.5f%c",a[i][j],' ');
      printf("\n");
   }
}

void PrintVec(a,title)
vector	a;
char	title[];
{
   int	i;

   printf("\n         %s \n",title);
   for (i=0;i<=2;i++)
      printf("%9.5f%c",a[i],' ');
   printf("\n");
}

void PrintVec10(a,title)
vector10  a;
char	  title[];
{
   int	i;

   printf("\n         %s \n",title);
   for (i=0;i<=4;i++)
       printf("%9.5f%c",a[i],' ');
   printf("\n");
   for (i=5;i<=9;i++)
       printf("%9.5f%c",a[i],' ');
   printf("\n");
   printf("\n");
}

void CopyMatrix (a,b)
matrix	a, b;
{
   short  i, j;

   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	 a[i][j] = b[i][j];
}

void MulMat(a,b,c)
matrix 	a, b, c;
{
   double  dummy;
   matrix  d, e;
   short i, j, k;

   CopyMatrix(d,b);
   CopyMatrix(e,c);
   for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
	 for (dummy=0, k=0; k<3; k++)
	    dummy += d[i][k] * e[k][j];
	 a[i][j] = dummy;
      }
}

void Transpose (a,b)
matrix	a, b;
{
   short i, j;
   matrix  c;

   CopyMatrix(c,b);
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	 a[j][i] = c[i][j];
}

void RotationMatrix(a,winkel,axis)
matrix 	a;
double	winkel;
char	axis;
{
   short  i, j, k;

   winkel = winkel * PI /180.0;
   for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	 a[i][j] = 0.0;
   switch (axis) {
   case 'x':
      i=0;
      j=1;
      k=2;
      break;
   case 'y':
      i=1;
      j=2;
      k=0;
      break;
   case 'z':
      i=2;
      j=0;
      k=1;
      break;
   }
   a[i][i] = 1.0;
   a[j][j] = cos(winkel);
   a[k][k] = cos(winkel);
   a[j][k] = sin(winkel);
   a[k][j] = -sin(winkel);
}

void NewLine(fp)
FILE  *fp;
{
   char  c;

   while ((c = getc(fp)) != '\n');
}


double Compliment(a)
double	a;
{
   if (a*a > 1.0) {
      printf(" Warning: Function \"Compliment\" called with a value greater than 1.0!! Returns a value of 0.0 ");
      return 0.0;
   }
   else
      return sqrt(1.0 - a*a);
}

void WriteSpec(low, high, sp, n,name)
double    low, high, *sp;
int	  n;
char	  name[];
{
   int    i;
   FILE   *fp, *fopen();
   double range;

   range = high - low;
   strncat(name,asc,sizeof(asc));
   if ((fp = fopen(name,"w")) == NULL)
      printf("file %s cannot be opened",name);
   for (i=0; i<n; i++)
      fprintf(fp,"%.3f\t%.9f\n",low + i*range/(n-1), sp[i]);
   fclose(fp);
}

