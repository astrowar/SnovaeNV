// snova.cpp : Defines the entry point for the console application.
//



#include <cstdio>
#include <cstdlib>
#include <cmath>


#include <vector>
#include <ctime>
#include <algorithm>

#define NP 400

#define Q 1.29
#define me 0.511
#define eta 0.0
#define Mk  2.14     //  Kamiokande detector mass (kton).
#define Mimb 6.8     // IMB detector mass (kton).
#define Mbaksan 0.28 // Baksan detector mass (kton).
#define fk 1.0
#define fimb 0.9055
#define fb 1.0
#define Cn    1/(4 * 3.1416)
#define lnVk   14.56
#define lnVimb  15.73
#define delt 1.0//2.3
#define timeK 20.0//10.43
#define timeIMB 16.0 //5.9






#define real double


using namespace std;

#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ; x+=0  ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };


#define Dump(  ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ;x+=(dx)   ){ printf(">%f %f \n",x, ( ff ) );};

#define H  6.0

//regiao de calculo da massa dos distibuicoes T
#define Tmin 1.5
#define Tmax 6.0
#define ddT  0.1

#define sigma_tp 1.0

//  tempo
#define tpmin 0.1
#define tpmax 4.5
#define ddtp 0.1 //2.0


// alpha
#define Amin 10.5
#define Amax 22.0
#define ddA  0.5

//  eps
#define epsmin 1.0
#define epsmax 50.0 // 50.0
#define ddeps  1.0 // 5.0


// tau1

#define  tau1min 0.05
#define  tau1max 5.0 // 30.0
#define  ddtau1  0.05 // 3.0



// tau2

#define  tau2min 0.0
#define  tau2max 4.0
#define  ddtau2  4.9


//Scale parameter
#define apmin  0.01
#define apmax  0.4
#define ddap   0.4


// Noise's parameters

#define pk     1.0
#define pimb   1.0
#define pb     1.0
#define effK   1e-5
#define effIMB 1e-5


//double PARMS[5];
//(FILE*)FILES[5];
//int NUMPARAM=5;





double Meff;


//vector <double> vAlpha, vAlpha2 ,vTemp , vDeltaTp , vAp ,vG  ;  //vetores que armazenam os parametros

vector <double> vAlpha, vAlpha2, vTemp, vtp, vAp, vtau1, vtau2;  //vetores que armazenam os parametros



double *****L;

double ******Lk; //vetor que armazena a likehood  K base para todos

double ****Limb; // vetor que armazena a likelihood IMB base para todos

double ******L_s12; //vetor que armazena a likehood S1 e S2

double *****L_g; //vetor que armazena a likehood g

double ***L_a; //vetor que armazena a likehood a

double Lmax;
//********************************************************** Data ***********************************************************************


//Kamikande events:


double tk[17] = { 0.0, 0.0, 0.107, 0.303, 0.324, 0.507, 0.686, 1.541, 1.728, 1.915, 9.219, 10.433, 12.439, 17.641, 20.257, 21.355, 23.814 }; // times of events;
double Ek[17] = { 0.0, 20.0, 13.5, 7.5, 9.2, 12.8, 6.3, 35.4, 21, 19.8, 8.6, 13, 8.9, 6.5, 5.4, 4.6, 6.5 };   // energy of events;
double Sigmak[17] = { 0.0, 2.9, 3.2, 2.0, 2.7, 2.9, 1.7, 8, 4.2, 3.2, 2.7, 2.6, 1.9, 1.6, 1.4, 1.3, 1.6 }; // standard deviation by events;
double Bk[17] = { 1, 1.6e-5, 1.9e-3, 2.9e-2, 1.2e-2, 2.1e-3, 3.7e-2, 4.5e-5, 8.2e-5, 1.5e-5,            // detector's noise;
1.5e-2, 1.9e-3, 1.6e-2, 3.8e-2, 2.9e-2, 2.8e-2, 3.8e-2 };


// IMB events:

double timb[9] = { 0.0, 0.0, 0.412, 0.650, 1.141, 1.562, 2.684, 5.010, 5.582 };
double Eimb[9] = { 0,38,37,28,39,36,36,19,22 }; // energy of events;
double Sigmaimb[9] = { 0,7,7 ,6,7,9, 6,5,5 };      // standard deviation by events;
double Bimb[9] = { 0,0,0,0,0,0,0,0,0 };         // detector's noise;


												// Defined functions and matrixes:



double**  gen_matrix2(int n, int m)
{
	double **d;

	d = (double**)malloc(sizeof(double*)*n);

	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = (double*)malloc(sizeof(double)*m);
	}

	return d;
}
double*** gen_matrix3(int n, int m, int m2)
{
	double ***d;

	d = (double***)malloc(sizeof(double**)*n);


	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = gen_matrix2(m, m2);
	}

	return d;
}

double**** gen_matrix4(int n, int m, int m2, int m3)
{
	double ****d;

	d = (double****)malloc(sizeof(double***)*n);

	for (int i3 = 0; i3<n; i3++)
	{
		d[i3] = gen_matrix3(m, m2, m3);
	}

	return d;
}

double***** gen_matrix5(int n, int m, int m2, int m3, int m4)
{
	double *****d;

	d = (double*****)malloc(sizeof(double****)*n);

	for (int i4 = 0; i4<n; i4++)
	{
		d[i4] = gen_matrix4(m, m2, m3, m4);
	}

	return d;
}

double****** gen_matrix6(int n, int m, int m2, int m3, int m4, int m5)
{
	double ******d;

	d = (double******)malloc(sizeof(double*****)*n);

	for (int i5 = 0; i5<n; i5++)
	{
		d[i5] = gen_matrix5(m, m2, m3, m4, m5);
	}

	return d;
}

double******* gen_matrix7(int n, int m, int m2, int m3, int m4, int m5, int m6)
{
	double *******d;

	d = (double*******)malloc(sizeof(double******)*n);

	for (int i6 = 0; i6 < n; i6++)
	{
		d[i6] = gen_matrix6(m, m2, m3, m4, m5, m6);
	}

	return d;
}



double  Integra_1(double *y, vector<double> x1)
{
	double soma = 0.0;
	int i;
	int ifinal = x1.size();
	for (i = 0; i< ifinal - 2; i += 2)
	{
		soma += (y[i] + 4 * y[i + 1] + y[i + 2])*(x1[i + 2] - x1[i]);
	}
	return soma / 6.0;
}


double  Integra_2(double **yf, vector<double> x1, vector<double> x2)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_1(yf[0], x2);
	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_1(yf[i + 1], x2);
		yc = Integra_1(yf[i + 2], x2);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;
}

double  Integra_3(double ***y, vector<double> x1, vector<double> x2, vector<double> x3)
{
	double soma = 0.0;
	double  ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_2(y[0], x2, x3);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_2(y[i + 1], x2, x3);
		yc = Integra_2(y[i + 2], x2, x3);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;


}

double  Integra_4(double ****y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_3(y[0], x2, x3, x4);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_3(y[i + 1], x2, x3, x4);
		yc = Integra_3(y[i + 2], x2, x3, x4);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;



}

double  Integra_5(double *****y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4, vector<double> x5)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_4(y[0], x2, x3, x4, x5);

	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_4(y[i + 1], x2, x3, x4, x5);
		yc = Integra_4(y[i + 2], x2, x3, x4, x5);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;



}


double  Integra_6(double ******y, vector<double> x1, vector<double> x2, vector<double> x3, vector<double> x4, vector<double> x5, vector<double> x6)
{
	double soma = 0.0;
	double   ya, yb, yc;
	int i;
	int ifinal = x1.size();
	ya = Integra_5(y[0], x2, x3, x4, x5, x6);
	for (i = 0; i< ifinal - 2; i += 2)
	{
		yb = Integra_5(y[i + 1], x2, x3, x4, x5, x6);
		yc = Integra_5(y[i + 2], x2, x3, x4, x5, x6);
		soma += (ya + 4 * yb + yc)*(x1[i + 2] - x1[i]);
		ya = yc;
	}
	return soma / 6.0;


}

 



double obs_sample_statistic_err(double mu)
{
	double soma;
	soma = 0;


	return soma;
}

double obs_sample_statistic(double mu, double dm)
{
	double soma;

	return soma;
}

//*************************************************** Defined functions ******************************************************************
  
 

 


typedef struct LikelihoodParameter
{
	real alpha, T, ap, tp, tau1, tau2;
	real result;

	LikelihoodParameter(real _alpha, real _T, real _ap, real _tp, real _tau1, real _tau2);
} LikelihoodParameter;

LikelihoodParameter::LikelihoodParameter(real _alpha, real _T, real _ap, real _tp, real _tau1, real _tau2)
{

	alpha = _alpha;
	T = _T;
	ap = _ap;
	tp = _tp;
	tau1 = _tau1;
	tau2 = _tau2;
	result = 0.0;
}




//**************************************************** Priori functions  ****************************************************************


double priori1(double alpha, double T, double tp, double ap, double tau1, double tau2) {




	//	return  1.0/ pow(alpha,2)  ;
	double pAlpha = 1.0 / pow(alpha, 2.0);
	double pTemp = 1.0 / pow(T, 4.0);
	//double pTau = 1.0 / (pow(tau1, 0.2) * pow(tau2, 0.2));


	return pTemp * pAlpha;


}



//**************************************************** Likelihood's Calculus *************************************************************


//**************************************************  Predictive's calculus  **********************************************************


// Prior 1 --> two gaussians:


double predictive1(void) {


	double soma = 1.0;



	return soma;






}

// Prior 2 --> 1 / M:

double predictive2(void) {


	return 0.0;




}

//********************************************** Probabilities calculus **************************************************************

//Probability 1:

double Probability1(int i1, int i2, int i3, int i4, int i5, int i6, int i7) {

	return 1;

}

void computeParams(std::vector<LikelihoodParameter> &params);
 



int build_Program()
{
	return  0;
};

int main() {



	FILE *arquivo;
	double i, x, y, z;

	double x1, x2, x3, x4, x5;
	double lmax1, lmax2, lmax3, lzmax;
	double predi;
	double ei1, ei2, ei3, ei4, ei5, ei6, ei7;

	double tmp_max = 0;

	printf("Start ! \n");

	Lmax = 0.0;


	ei1 = 0;
	ei2 = 0;
	ei3 = 0;
	ei4 = 0;
	ei5 = 0;
	ei6 = 0;
	ei7 = 0;

	srand(time(NULL));

	for (x = Amin; x <= Amax; x = x + ddA) vAlpha.push_back(x);

	for (x = Tmin; x <= Tmax; x = x + ddT) vTemp.push_back(x);


	for (x = tpmin; x <= tpmax; x = x + ddtp) vtp.push_back(x);

	for (x = apmin; x <= apmax; x = x + ddap) vAp.push_back(x);

	for (x = tau1min; x <= tau1max; x = x + ddtau1) vtau1.push_back(x);

	for (x = tau2min; x <= tau2max; x = x + ddtau2) vtau2.push_back(x);



	printf("Pts =%i  %i %i %i %i %i  \n", vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());
	printf("Size =%i \n", vAlpha.size()* vTemp.size()* vtp.size() *vAp.size() *vtau1.size() * vtau2.size());

	Lk = gen_matrix6(vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());
	//Limb = gen_matrix4(  vAlpha.size(), vTemp.size(), vDeltaTp.size() ,vAp.size() );


	//L_s12 = gen_matrix6(vAlpha.size(), vTemp.size(), vAp.size(), vtp.size(), vtau1.size(), vtau2.size());

	//L_s12= gen_matrix5( vAp.size(), vAlpha.size(), vTemp.size(), vtau1.size(),vtau2.size() );


	time_t t_inicial;
	time_t t_atual;
	t_inicial = time(NULL);


	std::vector<LikelihoodParameter> params;
	std::vector<LikelihoodParameter> paramsBuffer;
	LikelihoodParameter max_likelihood_parameter(0, 0, 0, 0, 0, 0);
	for (int i1 = 0; i1 < vAlpha.size(); i1++)
	{

		for (int i2 = 0; i2 < vTemp.size(); i2++)	 
			for (int i3 = 0; i3 < vAp.size(); i3++)
			{ 
				params.clear(); //Mantem um cache de memoria do ciclo anterior
				paramsBuffer.clear();
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)						
						{
						
							//printf(" %i %i %i %i %i \n", i1,i2,i3,i4,i5);


							//Lk[i1][i2][i3][i4][i5][i6] = Likelihood_combined(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							//printf("%6.2f %6.2f %6.2f %6.2f  %6.2f %6.2f : %15.12e \n", vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6], Lk[i1][i2][i3][i4][i5][i6]);

							//                 _alpha,      _T,       _ap,       _tp,     _tau1,      _tau2
							paramsBuffer.emplace_back(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							if (paramsBuffer.size() >= 1024 *2  )
							{
								computeParams(paramsBuffer); params.insert(params.end(), paramsBuffer.begin(), paramsBuffer.end()); paramsBuffer.clear(); 
								 
								 
							}
						 
						}
			 
				computeParams(paramsBuffer); params.insert(params.end(), paramsBuffer.begin(), paramsBuffer.end()); paramsBuffer.clear();
				 
				

				int iParam = 0;
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
						{
							double i_priori1 = priori1(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
							params[iParam].result = i_priori1 * params[iParam].result;
							Lk[i1][i2][i3][i4][i5][i6] = params[iParam].result;

							if (params[iParam].result > Lmax)
							{
								Lmax = params[iParam].result;
								max_likelihood_parameter = params[iParam];
							}
							iParam++;

							// if((iParam%300) ==0 ) printf("%6.2f %6.2f %6.2f %6.2f  %6.2f %6.2f : %15.12e \n", vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6], Lk[i1][i2][i3][i4][i5][i6]);
						}

			}

		t_atual = time(NULL);
		printf(" %i  %i   restam %4.1f  minutos \n", i1, vAlpha.size(), (vAlpha.size() - (i1 + 1))*(((t_atual - t_inicial) / 60.0) / (i1 + 1)));

	}



	predi = Integra_6(Lk, vAlpha, vTemp, vAp, vtp, vtau1, vtau2);


	// Normaliza a preditivTemp, aplica a priori, e translada as matrizes


	lzmax = 0.0;
	for (int i1 = 0; i1 < vAlpha.size(); i1++) {
		for (int i2 = 0; i2 < vTemp.size(); i2++)
			for (int i3 = 0; i3 < vAp.size(); i3++)

			{
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
							break;

				//Lk[i1][i2][i3][i4][i5][i6] = Lk[i1][i2][i3][i4][i5][i6] /  predi;

				//Lk[i1][i2][i3][i4][i5][i6] = Lk[i1][i2][i3][i4][i5][i6] * priori1(vAlpha[i1], vTemp[i2], vAp[i3], vtp[i4], vtau1[i5], vtau2[i6]);
				//	Limb[i1][i3][i4][i5] = Limb[i1][i3][i4][i5] * priori1(vAlpha[i1],vTemp[i3],vDeltaTp[i4], vAp[i5] );

			}
	}

	lzmax = 0.0;
	lmax1 = vAlpha[0];

	lmax3 = vtp[0];

	FILE *fmm;


	fmm = fopen("parameters.dat", "w+");
	fprintf(fmm, "Max Likehood Parameter:%g\n alpha %6.2f\n T %6.2f\n Tp %6.2f\nAmpliture %6.2f\n Tau 1 %6.2f\n Tau 2%6.2f \n ", max_likelihood_parameter.result,
		max_likelihood_parameter.alpha, max_likelihood_parameter.T,
		max_likelihood_parameter.tp,
		max_likelihood_parameter.ap,
		max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);

	fclose(fmm);

	printf("Finalizado\n ");


	fmm = fopen("Lat.dat", "w+");


	for (int i1 = 0; i1 < vAlpha.size(); i1++)
	{
		for (int i2 = 0; i2 < vTemp.size(); i2++)

		{
			double soma = 0;


			for (int i3 = 0; i3 < vAp.size(); i3++)
				for (int i4 = 0; i4 < vtp.size(); i4++)
					for (int i5 = 0; i5 < vtau1.size(); i5++)
						for (int i6 = 0; i6 < vtau2.size(); i6++)
							soma += ddap * ddtp * ddtau1 * ddtau2 * Lk[i1][i2][i3][i4][i5][i6];

			fprintf(fmm, "%f %f %g \n", vAlpha[i1], vTemp[i2], soma);

		}

		fprintf(fmm, "\n");
	}

	fclose(fmm);



	fmm = fopen("tau12.dat", "w+");


	for (int i5 = 0; i5 < vtau1.size(); i5++)
	{
		for (int i6 = 0; i6 < vtau2.size(); i6++)

		{
			double soma = 0;
			for (int i1 = 0; i1 < vAlpha.size(); i1++)
				for (int i2 = 0; i2 < vTemp.size(); i2++)
					for (int i3 = 0; i3 < vAp.size(); i3++)
						for (int i4 = 0; i4 < vtp.size(); i4++)

						{
							soma += ddap * ddtp * ddA* ddT*   Lk[i1][i2][i3][i4][i5][i6];
						}
			fprintf(fmm, "%f %f %g \n", vtau1[i5], vtau2[i6], soma);
		}

		fprintf(fmm, "\n");
	}

	fclose(fmm);







	fmm = fopen("Lap.dat", "w+");



	for (int i4 = 0; i4 < vtp.size(); i4++)
	{
		for (int i3 = 0; i3 < vAp.size(); i3++)
		{
			double soma = 0;
			for (int i6 = 0; i6 < vtau2.size(); i6++)
				for (int i5 = 0; i5 < vtau1.size(); i5++)
					for (int i2 = 0; i2 < vTemp.size(); i2++)
						for (int i1 = 0; i1 < vAlpha.size(); i1++)
						{
							soma += ddtau1 * ddtau2*ddT * ddA * Lk[i1][i2][i3][i4][i5][i6];
						}

			fprintf(fmm, "%f %f %g \n", vtp[i4], vAp[i3], soma);

		}
		fprintf(fmm, "\n");
	}




	fclose(fmm);


	printf("Likelihood = %g \n", Lmax);

	fmm = fopen("rcol.dat", "w+");
	for (double e = 1; e < 50; e += 0.2)
	{
		for (double Te = 0.1; Te < 20; Te += 0.2)
		{
			//fprintf(fmm, " %f %f %e \n", e, Te, etabarK(e)*Rcol(e, max_likelihood_parameter.alpha, 0, Te, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2));
		}
		fprintf(fmm, "\n");
	}
	fclose(fmm);


	fmm = fopen("temp.dat", "w+");
	{
		for (double tt = 0; tt < 20; tt += 0.02)
		{
			//double Te = Temp(tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
			//fprintf(fmm, " %f %f \n", tt, Te);
		}

	}
	fclose(fmm);


	fmm = fopen("alpha.dat", "w+");
	{
		for (double tt = 0; tt < 20; tt += 0.02)
		{
			//double Te = get_alpha(max_likelihood_parameter.alpha, tt, max_likelihood_parameter.T, max_likelihood_parameter.ap, max_likelihood_parameter.tp, max_likelihood_parameter.tau1, max_likelihood_parameter.tau2);
			//fprintf(fmm, " %f %f \n", tt, Te);
		}

	}
	fclose(fmm);

	return 0;

}






