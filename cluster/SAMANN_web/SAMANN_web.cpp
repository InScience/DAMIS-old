
// samann_mpi_nuos.cpp : Defines the entry point for the console application.
//prdiniai svoriai imami ish failu: w ir ww
//universali programa test duom aibems
//failo pavadinimas test_duom.txt
//
//#include "stdafx.h"
#include <iomanip>
#include <iostream.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <fstream.h>
#include <mpi.h>
#include <string>
using namespace std;
#include "tools.h"
#include <mpi.h>

char ats;
int K, M, L;

long	_rows, _cols;
double** _pData;

double d1[6][410][410], alpha, alpha_tarp;
double d2[6][410][410];
double delta_w[410][410];
double delta_ww[410][6];
double b[200000][200], aa[200000][200], a[200000][200];
double y[200000][200];
double z[200000][11], v_senas[100][100], vv_senas[100][10];
double momentum;
long   n,n2, n12, i, j, d, ii, ii2, ii_new, jj, nn1, nn2, k,itt;
long   g, fil, pr_sk, IT;
double _max[200];
double _min[200] ,laikas;
int    r[100000], mas1[100000];
double min_klaida;

//double v[210][210];
double v[100][100], v_temp[3000], vv_temp[2000], b_temp[800000];
//double vv[210][210];
double vv[100][10];
double s,x, eta, eta1, eta2, l, s_kl2, ss1, ss2, sam_min0, sam_min1, sam_min2;
double delta1[500][500];
double delta2[500], sum[500];
int  it_sk1, it_sk2;
double s_kl[3];

long   sl, w, p, m, ww;
long   v1, v2;
/*d - vektoriaus elementu sk;
     nn1,nn2 - du vektoriai;
     n - vektoriu skaicius; sl - pasl. sluoksnio neuronu sk.; */


time_t start1, start2, finish1, finish2, finish3;

FILE *file1, *file6, *file7;
//const char *f1 = "data.txt";

const char *f6 = "sv01.txt";
const char *f7 = "svv01.txt";

char *pFileName = new char[100];
   char *pResultsFile = new char[100];


void sammon_l(int);
void sammon2(int, int, double, int, double, double);
void sammon3(int, int, double, int, double, double);
void sammon(int nn, double alpha);
void sammon_b(int nn, double alpha);
double paklaida(void);

void sammon_l(int n)
{
  int  i, j, k;
  double  s;
  sum[1]=0;
  s=0;

  for (i=1; i<=n-1; i++)
    {
      for (j=i+1; j<=n; j++)
	  {
		sum[2]=0;
		for (k=1; k<=d; k++)
		{
			s=(a[i][k]-a[j][k])*(a[i][k]-a[j][k]);
			sum[2]=sum[2]+s;
			
		}
		//printf("tarp atst:  %10.7lf \n", sum[2]);
		sum[1]=sum[1]+sqrt(sum[2]);
	  }	
    }
    l=1/sum[1];
}


void sammon2(int nn1, int nn2, double momentum, int ww, double eta, double alpha)
{
  int i, j, k;
  double tr, s, d_proj, d_pr;
  double sum[5];
  sum[1]=0;
  sum[3]=0;
  for (i=1; i<=d; i++)
    {
       s=(a[nn1][i]-a[nn2][i])*(a[nn1][i]-a[nn2][i]);
       sum[1]=sum[1]+s;	
	}
  d_pr=sqrt(sum[1]);

  if (d_pr==0)  d_pr=0.00001;

  sum[2]=0;
  for (k=1; k<=2;k++)
	{
      s=(z[nn1][k]-z[nn2][k])*(z[nn1][k]-z[nn2][k]);
      sum[2]=sum[2]+s;
    }
	d_proj=sqrt(sum[2]);

    if (d_proj==0) d_proj=0.00001;


  tr=(d_pr-d_proj)/(d_pr*d_proj);
  delta2[1]=-2*l*tr*(z[nn1][1]-z[nn2][1]); 
  delta2[2]=-2*l*tr*(z[nn1][2]-z[nn2][2]); 

    
  for (j=0; j<=sl; j++)
    {
      for (k=1; k<=2; k++)
		{                
			d2[1][j][k]=delta2[k]*(1-z[nn1][k])*z[nn1][k];
			d2[2][j][k]=delta2[k]*(1-z[nn2][k])*z[nn2][k];
			if (ww==1) vv_senas[j][k]=0;
				else 	vv_senas[j][k]=delta_ww[j][k];
			delta_ww[j][k]=-eta*alpha*(d2[1][j][k]*y[nn1][j]-d2[2][j][k]*y[nn2][j]);		
			vv[j][k]=vv[j][k]+delta_ww[j][k]+momentum*vv_senas[j][k];
	  }  
     }
 }


void sammon3(int nn1, int nn2, double momentum, int ww, double eta, double alpha)
{
  int i,j, k;
  double x, xx;
  double s[5];
  x=0;
  xx=0;
  for (j=0;j<=sl; j++)
    {
      s[1]=0;
      s[2]=0;
      for (k=1; k<=2; k++)
		{
          x=d2[1][j][k]*vv[j][k];
          s[1]=s[1]+x;
          xx=d2[2][j][k]*vv[j][k];
          s[2]=s[2]+xx;
		}
      delta1[1][j] = s[1];
      delta1[2][j] = s[2];
    }

  for (i=0; i<=d; i++)
    {
      for (j=1; j<=sl; j++)
		{
			d1[1][i][j]=delta1[1][j]*(1-y[nn1][j])*y[nn1][j];
		    d1[2][i][j]=delta1[2][j]*(1-y[nn2][j])*y[nn2][j];
			if (ww==1) v_senas[i][j]=0;
				else v_senas[i][j]=delta_w[i][j];
			delta_w[i][j]=-eta*alpha*(d1[1][i][j]*a[nn1][i]-d1[2][i][j]*a[nn2][i]);
			v[i][j]=v[i][j]+delta_w[i][j]+momentum*v_senas[i][j];
			
		}
    }
}


double paklaida(void)
{
	double sum[9], tk[3];
	double err, ss_a, ss_z;
	int k, m;
	

		sum[5]=0;
		sum[6]=0;
		sum[3]=0;
		sum[4]=0;
        for (i=1; i<=(n-1); i++)
		{
		sum[3]=0;
		sum[4]=0;
	    for (j=(i+1); j<=n; j++)
			{
				sum[2]=0;
				sum[1]=0;
			    for (k=1; k<=d; k++)
					{
						ss_a=(b[i][k]-b[j][k])*(b[i][k]-b[j][k]);
						sum[1]=sum[1]+ss_a;
					}
				if (sum[1]==0) sum[1]=0.00001;
				//vektoriai sutampa. koreguojam kintamuosius, kad nebutu dalybos is nulio
				tk[1]=sqrt(sum[1]);
        			sum[3] = sum[3] + tk[1];
	
				for (m=1; m<=2; m++)
					{
					ss_z=(z[i][m]-z[j][m])*(z[i][m]-z[j][m]);
					sum[2]=sum[2]+ss_z;
					}
				tk[2]=sqrt(sum[2]);
				sum[4] = sum[4] + ((tk[1]-tk[2])*(tk[1]-tk[2]))/tk[1];
		//printf(" sum1 sum2  sum4 %10.7lf %10.7lf %10.7lf \n", tk[1], tk[2],  ((tk[1]-tk[2])*(tk[1]-tk[2]))/tk[1]);


//printf(" error  %10.7lf \n", (tk[1]-tk[2]));
			}

			sum[5]=sum[5]+sum[3];
			sum[6]=sum[6]+sum[4];
		}
    		err = (1/sum[5])*sum[6];

				
			return err;
}


void sammon(int nn,double alpha)
{
  double xx, x;
  int    i,j,k;
  double sum[6];
  for (j=1; j<=sl; j++)
    {
      sum[1]=0;
      for (i=0; i<=d; i++)
	  { 
		  x=v[i][j]*a[nn][i];
		  sum[1]=sum[1]+x;
	  }
	   y[nn][j] = 1/(1+(exp(alpha*(-sum[1]))));
    }

  for (k=1; k<=2; k++)
    {
      sum[2]=0;
      for (j=0; j<=sl; j++)
		{
		  xx=vv[j][k]*y[nn][j];
		  sum[2]=sum[2]+xx;
		}

      z[nn][k] = 1/(1+(exp(alpha*(-sum[2]))));

    }
}


void sammon_b(int nn, double alpha)
{
  double xx, x;
  int    i,j,k;
  double sum[6];
  for (j=1; j<=sl; j++)
    {

      sum[1]=0;
      for (i=0; i<=d; i++)
	  { 
		  x=v[i][j]*b[nn][i];
		  sum[1]=sum[1]+x;
	  }
	   y[nn][j] = 1/(1+(exp(alpha*(-sum[1]))));

    }

  for (k=1; k<=2; k++)
    {
      sum[2]=0;
      for (j=0; j<=sl; j++)
		{
		  xx=vv[j][k]*y[nn][j];
		  sum[2]=sum[2]+xx;
		}
      z[nn][k] = 1/(1+(exp(alpha*(-sum[2]))));

    }
}

int main(int argc, char **argv)
{
 int K,L;
 int pp;
  int myid;
  MPI_Status stat;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &pp);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

if (myid>=0) 
  
{

it_sk2=atoi(argv[2]);   //iteraciju skaicius
//eta1=atoi(argv[3]);   //iteraciju skaicius
pFileName=argv[1];
pResultsFile=argv[4];


//d= 4;  /*vektoriaus elementu skaicius*/
//n= 150; /*vektoriu skaicius*/
//n2=100; /*duomenu aibes poaibis*/
sl= 10; /*pasleptojo sluoksnio el. skaicius  sl*/
eta1=10;
it_sk1=1000;
//it_sk2=1000;

	  s_kl[1]=0;
      sam_min1=1000;
}


if (myid  >= 0)
    {

      file6 = fopen(f6, "r");
      file7 = fopen(f7, "r");

   double num;
   ifstream file;
 
file.open(pFileName, ifstream::in);

    if (file.is_open()) 
	{
        //Nunuliname pradinius masyvo matmenis
        _rows = -1;
        _cols = 1000;
        // Pasiimame masyvo dydi
        j = 0;
        while (!file.eof()) 
		{	
            i = 0;
            string line;
            getline(file, line, '\n');
            char *tmpstr, *tmpstr2;
            tmpstr = &line[0];
            while (strlen(tmpstr) > 0 && i<_cols) 
			{
                num = double (strtod(tmpstr, &tmpstr2));
                tmpstr = tmpstr2;
                i++;
            }

            // Pirma gerai nuskaityta eilute jei sulpeliu maziau negu 1 arba nesutampas su cols skaiciumi
            if (j == 0) 
			{
                if (i > 1) { _cols = i; } else { j = -1; }
            } else 
			{
                if (i != _cols) { j--; }
            }
            j++;
        }
        _rows = j;

	}  


d=_cols;
n=_rows;

//n2=long(n/3);
//n2=100;

//if (n>=50000) n2=long(n/50); else if ((n>=10000)&&(n<40000)) n2=long(n/40); else if ((n>=3000)&&(n<1000)) n2=long(n/30); else if (n>=3000) n2=long(n/15); else if ((n>1000)&&(n<3000)) n2=long(n/10); else if ((n>=300)&&(n<=1000)) n2=long(n/5); else n2=long(n/3);
if (n>=50000) n2=500; else if ((n>=10000)&&(n<50000)) n2=300; else if ((n>=1000)&&(n<10000)) n2=300; else n2=long(n*0.5);


const char *f1 = pFileName;
   file1 = fopen(f1, "r");   


/*
for (i=0; i<n; i++)
				{

			      for (j=0; j<d	; j++)
					{   

					  	fscanf(file1, "%lf", &b[i+1][j+1]);  //file1
					}                                 
					fscanf(file1, "\n");  //file1
				}	

*/


// laikinai

     if (_rows > 0 && _cols > 0) 
	{
        // pakartotinai nuskaitome faila dabar jau zinodami masyvo dydi
        // isskiriam vietos duomenims
    
		 //  _pData = (double*) realloc (_pData, _rows * _cols * sizeof(double));
	  _pData=(double**) malloc(_rows* sizeof(double*));
	  for (i=0; i<=_rows; i++) 
	  {
		  _pData[i]=(double*) malloc(_cols* sizeof(double));
	  }

//         file.open("test_duom_iris.txt", ifstream::in);
        if (file.is_open()) 
		{

            j = 0;
            while (!file.eof() && j < _rows) 
			{
                i = 0;
                string line;
                getline(file, line, '\n');
                char *tmpstr, *tmpstr2;
                tmpstr = &line[0];
                while (strlen(tmpstr) > 0 && i < _cols) 
				{
                    num = double(strtod(tmpstr, &tmpstr2));
                    _pData[i][j]=num;
			printf("_pData     %12.8lf \n", _pData[i][j]);
                     printf("test test test test \n");
			//*(_pData + j*_cols + i) = num;
                    tmpstr = tmpstr2;
                    i++;
                }

                // Pirma gerai nuskaityta eilute jei sulpeliu maziau negu 1 arba nesutampas su cols skaiciumi
                if (j == 0) 
				{
						if (i < _cols) 
						{
						  j = -1;
						}
                } else 
				{
                    if (i != _cols) 
					{
                        j--;
                    }
                }
                j++;
            }
  
        }
	 }
  





for (i=0; i<n; i++)
				{
			      for (j=0; j<d	; j++)
					{   
					  	b[i+1][j+1]=_pData[i][j];
					}                                 
				}	

// */ 


file.close();
srand( (unsigned)time( NULL ));
	 
	if (file6 == NULL) printf("Failas neatidarytas");
      else
		{
			for (i=0; i<=d; i++)
				{
			      for (j=1; j<=sl; j++)
					{   
					fscanf(file6, "%lf", &v[i][j]); 

					}                                 
					fscanf(file6, "\n");

				}	
		}


	  
	  if (file7 == NULL) printf("Failas neatidarytas");
      else
		{
			for (i=0; i<=sl; i++)
				{
			      for (j=1; j<=2; j++)
					{  
					  fscanf(file7, "%lf", &vv[i][j]); 
				  }                                 
					  fscanf(file7, "\n");
				}	
		}






      _max[1]=0;
      _min[1]=10000;
      sum[1]=0;
      for (i=1; i<=(n-1); i++) 
		{
			for (j=(i+1); j<=n; j++)
				{
					sum[1]=0;
					for (k=1; k<=d; k++)
						{
							s=(b[i][k]-b[j][k])*(b[i][k]-b[j][k]);
							sum[1]=sum[1]+s;
						}
					if (sqrt(sum[1])>_max[1])   _max[1]=sqrt(sum[1]);
					if (sqrt(sum[1])<_min[1])   _min[1]=sqrt(sum[1]);
				}
		}

      for (j=1; j<=n; j++)  /*normuojam ish karto visus vektorius*/
		{ 
			for (k=1; k<=d; k++)
				{
					b[j][k]=b[j][k]/_max[1]; //((_max[1]-_min[1])*sqrt(d/2));
                 }
		
		}
}


 alpha=5;
 
IT=it_sk2;
time(&start1);
iteracija:


  long ats_sk;
 for (i=1; i<=n; i++)              
    mas1[i]=i;
 	L=n; 	
	M=0;
	srand( (unsigned)time( NULL )+myid);
	
gen2:
	ats_sk=1+(rand()%L);
	M=M+1;	
	r[M]=mas1[ats_sk];	
	for (i=ats_sk; i<= (L-1); i++)	
	mas1[i]=mas1[i+1];
	L=L-1;
	if (M<n2) goto gen2;

	for (i=1; i<=n2; i++)
		for (j=1; j<=d; j++)
{
		a[i][j]=b[r[i]][j];
}



      for (i=1; i<=n; i++) 
		{

			b[i][0]=1;
			y[i][0]=1;
			z[i][0]=1;
		}
 for (i=1; i<=n2; i++) 
		{
			a[i][0]=1;

		}

  sammon_l(n2);

  
  //Skaiciuojam pirmos aibes A1 isejimo vektorius ir svorius 
 // for (w=1; w<=it_sk1; w++)
w=1;
iteracija2:
		{
if (w>=2) alpha=1;
ww=1;
		for (v1=1; v1<=(n2-1); v1++)
					{
					for (v2=(v1+1); v2 <= n2; v2++)
						    {
								sammon(v1, alpha);
								sammon(v2, alpha);
								sammon2(v1, v2, 0.05, ww, eta1, alpha);
								sammon3(v1, v2, 0.05, ww, eta1, alpha);
								ww=ww+1;

				   	  		}
				
					}

			//if (w%10==0)
			{
				for (i=1; i<=n; i++)
					sammon_b(i, alpha);
					//double s_kl[3];
					/* !!!!!!!!!!!  */
					s_kl[1]=paklaida();
					if (s_kl[1]<sam_min1) 
						{
							sam_min1=s_kl[1];
							itt=w;
						}	

					time(&finish1); //fiksuojam laika A1
					
					

				
			}
     //printf("procesorius: %8ld    %10.7lf   %10.7lf  Nr.:  %ld   %10.1f \n", myid, s_kl[1], sam_min1, w, difftime(finish1, start1));
	laikas=difftime(finish1, start1);






//	if (difftime(finish1, start1)<54000) {w=w+1; goto iteracija2;}    // kol kas laiko nereikia
//}	 
}	  


IT=IT-1;
if (IT>0) { goto iteracija;}

	
    sprintf(pResultsFile, "%s", argv[4]);
    char *pRank = new char[4];
    sprintf(pRank, "%d", myid);
    strcat(pResultsFile, pRank);
    strcat(pResultsFile, ".txt");
   // pMds->SaveResults(pResultsFile, time);


//ifstream file;
//file.open(pResultsFile, ifstream::in);

ofstream file;
    file.open (pResultsFile);
if (file.is_open()) {

      file << "[Method]:" << "SAMANN" << "\n";
      file << "[InitMethod]:" << "none" << "\n";
      file << "[Parameters]:" << "it. skaicius:" <<it_sk2<< "\n";
      file << "[DataSet]:" << d << "x" << n << "\n";
      file << "[T]:" << fixed << setprecision(3) << laikas << "\n";
      file << "[I]:" << it_sk2 << "\n";
      file << "[E]:" << fixed << setprecision(15) << s_kl[1] << "\n";
      file << "[MWE]:" << fixed << setprecision(15) << "none" << "\n";
      file << "[SC]:" << fixed << setprecision(15) << "none" << "\n";

        for (int i=1; i<=n; i++) {
             file << z[i][1]<<";"<<z[i][2]<<";3"<<"\n";

        }
        file.close();
    } else {
        printf("negali atidaryti failo!\n");
    }
		
	 
		fclose(file1);
		fclose(file6);
		fclose(file7);
	
		MPI_Finalize();
		return 0;
}
