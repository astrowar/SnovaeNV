#include "cuda_runtime.h"

#include "device_launch_parameters.h"

#include <stdio.h>
#include <vector>
#define real double
#define number float


#define NP 400

#define Q 1.29
#define me 0.511
#define eta 0.0
#define Mk  2.14     //  Kamiokande detector mass (kton).
#define Mimb 6.8     // IMB detector mass (kton).
#define Mbaksan 0.28 // Baksan detector mass (kton).
#define fk (1.0)
#define fimb (0.9055)
#define fb (1.0)
#define Cn    (1.0/(4 * 3.1416))
#define lnVk   14.56
#define lnVimb  15.73
#define delt 1.0  //2.3
#define timeK 20.0  //10.43
#define timeIMB 16.0   //5.9
#define timeEnd 40.0   



#define epsmin 1.0
#define epsmax 50.0 // 50.0
#define ddeps   0.1 // 5.0


// Noise's parameters

#define pk     1.0
#define pimb   1.0
#define pb     1.0
#define effK   1e-5
#define effIMB 1e-5

#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ; x+=0  ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };
#define H  6.0

typedef struct LikelihoodParameter
{
	double alpha, T, ap, tp, tau1, tau2;
	double result;

	LikelihoodParameter(double _alpha, double _T, double _ap, double _tp, double _tau1, double _tau2);
} LikelihoodParameter;
typedef struct PressSchecter
{
	number ap, tp, tau1, tau2;
}PressSchecter;


__device__ number gaussian(number  x, number x0, number sigma)
{
	if (fabs(x - x0) < sigma * H) return (1.0 / (sqrt((number)2 * 3.1415928* sigma*sigma))) * exp(-0.5*pow((x - x0) / sigma, (number)2.0));
	return 0.0;
}



// Function type Press - Schecter
__device__  number Temp(number tu, number T, PressSchecter tmp) {


	if (tu > tmp.tp)
	{
		return  T * exp(-1.0 * min((((tu - tmp.tp) / tmp.tau1)), (number)12.0));
	}
	return  T;

}


__device__  number r(number T) {
	return 1.0;
}


__device__ number f_fermi(number E, number T)
{
	//if (T < 0.01) return 0;
	number ET = min(15.0,E / T);
	 
	return     1.0 / (exp(ET) + 1.0);
}


__device__ number kappa(number E) {

	const number a = 1.0 - Q / E;
	const number b = 1.0 - (2.0 * Q / E);
	const number c = (pow((number)Q, (number)2.0) - pow((number)me, (number)2.0)) / pow(E, (number)2.0);
	return a * sqrt(max(b + c, (number)0.0));
}


// Rate's neutrino - cooling component=

__device__  number Rcol(number E, number alpha, number T, number MMeff) {

	number fm;
	number kp;
	number saida = 0; 
	fm = f_fermi(E, T);
	if (fm <= 0.0) { return 0.0; }
	number alpha_t = alpha;
	kp = kappa(E);
	if (kp <= 0.0) return 0.0;
	saida = (1.22e-5) * pow(alpha_t, (number) 2.0) * MMeff * (pow(E, (number)4.0)) * fm * kp * pow(r(T), (number)2.0);
	 
	return saida;

}

__device__  number etabarK(number E) 
{
	 
	number c = 0.95*(1.0 - exp(-pow((E) / 9.3, 4.0)));
	return max(c,0.0);
	
}


__device__  number  noiseK(number E)
{
	return  effK * (gaussian(E, 6.0, 1.0) + 0.001); 
}

__device__  number StepK(number eps) {
	if (eps >  5.0) {
		return 1;
	}
	else { return 0.0; }
}

 

__device__  real LikelihoodK(number alpha, number T, PressSchecter tmp, real *LMax) {


	const number    tk[17] = { 0.0, 0.0, 0.107, 0.303, 0.324, 0.507, 0.686, 1.541, 1.728, 1.915, 9.219, 10.433, 12.439, 17.641, 20.257, 21.355, 23.814 }; // times of events;
	const number    Ek[17] = { 0.0, 20.0, 13.5, 7.5, 9.2, 12.8, 6.3, 35.4, 21, 19.8, 8.6, 13, 8.9, 6.5, 5.4, 4.6, 6.5 };   // energy of events;
	const number    Sigmak[17] = { 0.0, 2.9, 3.2, 2.0, 2.7, 2.9, 1.7, 8, 4.2, 3.2, 2.7, 2.6, 1.9, 1.6, 1.4, 1.3, 1.6 }; // standard deviation by events;
	const number    Bk[17] = { 1, 1.6e-5, 1.9e-3, 2.9e-2, 1.2e-2, 2.1e-3, 3.7e-2, 4.5e-5, 8.2e-5, 1.5e-5,            // detector's noise;
		1.5e-2, 1.9e-3, 1.6e-2, 3.8e-2, 2.9e-2, 2.8e-2, 3.8e-2 };
	 
	int  i ;
	number soma;


	real termo1, termo2;
	real prod;
	number eps;
	 

	number e1, e2;

	soma = 0.0;
	 
	 
	//number jddtp = 0.05;
	number time_end = timeEnd;
	termo1 = 0.0;
	time_end = min(3*timeEnd, tmp.tp + H * 1.0/tmp.tau1);
	number jddtp = (time_end - tmp.tp) / 30.0;
	//return (etabarK(eps)*(Cn*Rcol(eps + Q, alpha, Temp(0.0, T, tmp), Mk) + noiseK(epsmin)));
	
	{
		number Tj = Temp(0, T, tmp); 
		for (eps = epsmin; eps <= epsmax; eps = eps + ddeps)
		{
			number _etabarK = etabarK(eps); 
			termo1 += (_etabarK*(Cn*Rcol(eps + Q, alpha, Tj, Mk) + noiseK(eps)));
			 
		}
		termo1 = termo1 * tmp.tp;
	}

	for (number ti = tmp.tp; ti <= time_end; ti = ti + jddtp)
	{
		number Tj = Temp(ti, T, tmp);
		if (Tj > 0.01)
		{
			for (eps = epsmin; eps <= epsmax; eps = eps + ddeps)
			{
				termo1 += (etabarK(eps)*(Cn*Rcol(eps + Q, alpha, Tj, Mk) + noiseK(eps)))   * jddtp ;
			} 
			//Integra(termo1, (etabarK(eps)*(Cn*Rcol(eps + Q, alpha, Tj, Mk) + noiseK(eps))), eps, epsmin, epsmax, ddeps);
		}
	}
	termo1 = termo1 * ddeps;
	 
	//{
	//	number ti = 0;
	//	Integra(termo1, ProbNotIntegrate(alpha, ti, T, tmp), ti, 0, time_end, jddtp);
	//}

	prod = 1.0;
	for (i = 1; i <= 12; i++)
	{
		if (i == 6) continue;
		termo2 = 0.0;

		number Tj = Temp(tk[i], T, tmp);
		if (Tj > 0.01)
		{
			e1 = max(epsmin, Ek[i] - H * Sigmak[i]);
		    e2 = min(epsmax, Ek[i] + H * Sigmak[i]);
		    number dSigma = Sigmak[i] / 4.0;
			//Integra(termo2, RDet(eps) , eps, e1, e2, dSigma);
			Integra(termo2, (StepK(eps)* Cn * lnVk * gaussian(eps, Ek[i], Sigmak[i])*Cn*(Rcol(eps + Q, alpha, Tj, Mk) + 1.0*noiseK(eps))), eps, e1, e2, dSigma);
		}
		 
		prod = prod * (Bk[i] + termo2);
	}

	soma = exp(-termo1) *  prod;

	 
	if (soma > (*LMax)) *LMax = soma;
	return soma;
	
}





__device__ real Likelihood_combined(real alpha, real T, real ap, real tp, real tau1, real tau2)
{

	PressSchecter psch = { (number)ap, (number)tp, (number)tau1,(number)tau2 };
	real LMax = 0.0;
	return LikelihoodK(alpha, T, psch, &LMax);

}
 
__global__ void Likelihood(const LikelihoodParameter *inputParams,   double *results , unsigned int size)
{
	//int tid = threadIdx.x;
	unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < size)
	{
		results[tid]   = Likelihood_combined(inputParams[tid].alpha, inputParams[tid].T, inputParams[tid].ap, inputParams[tid].tp, inputParams[tid].tau1, inputParams[tid].tau2);
		 
	}
	__syncthreads();

}


bool IsPowerOfTwo(unsigned int x)
{
	return (x & (x - 1)) == 0;
}
// Helper function for using CUDA to add vectors in parallel.
cudaError_t LikelihoodList(LikelihoodParameter *a, double *results, const unsigned int size)
{
	cudaError_t cudaStatus;
	
	static LikelihoodParameter *dev_a = 0;
	static double *dev_Lk = 0;
 
	unsigned int  size_truc = size;

	int blockSize;      // The launch configurator returned block size 
	int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch 


	while (size_truc % 32 != 0) size_truc++;
	 
	int nth = size_truc / 64 + 1;

 

 
	if (dev_a == 0)
	{
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		cudaStatus = cudaMalloc((void**)&dev_a, size_truc * sizeof(LikelihoodParameter));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc  A failed! %i \n", size_truc * sizeof(LikelihoodParameter));
			goto Error;
		}
		cudaStatus = cudaMalloc((void**)&dev_Lk, size_truc * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc B failed! %i \n", size_truc * sizeof(double));
			goto Error;
		}
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(LikelihoodParameter), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed! \n");
		goto Error;
	}

 
	       // The actual grid size needed, based on input size 
 
	cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, Likelihood, 0, size_truc);
	//  gridSize = (size_truc + blockSize - 1) / blockSize;

	//printf("gridSize = %i ,blockSize = %i \n", gridSize, blockSize);
	// Launch a kernel on the GPU with one thread for each element.
	  nth = size_truc / 64 + 1;

	//printf("BLOCK = %i ,THREADS_PER_BLOCK = %i \n", nth, 64);
	Likelihood <<<nth, 64 >>>(  dev_a  , dev_Lk , size );
	//Likelihood << <gridSize, blockSize >> >(dev_a, dev_Lk, size);

	cudaDeviceSynchronize();

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "size trunc: %u\n", size_truc);
		fprintf(stderr, "Likelihood launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(results, dev_Lk, size * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	 

Error:
	 
	//cudaFree(dev_a);
	//cudaFree(dev_Lk);

	cudaDeviceSynchronize();
	return cudaStatus;
}




void computeParams(std::vector<LikelihoodParameter> &params)
{
	static int enQueues = 0;
	unsigned int size = params.size();
	//printf("Queues to compute  = %i \n", size);

	std::vector<double> results(size, 0.0);	
	LikelihoodList(params.data(), results.data(), size);
	for(unsigned int k = 0;k< size ;++k)
	{
		//if (k==0) printf("result  = %g  \n" ,results[k]);
		params[k].result = results[k];
	}
	enQueues += size;

}