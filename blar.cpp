#include <vector>
#include <iostream>

extern "C"
{
	void dsaupd_(int*, const char*, const int*, const char*, const int*, const double*, double*, const int*, double*, const int*, int*, int*, double*, double*, const int*, int*);
	void dseupd_(const bool*, const char*, bool*, double*, double*, const int*, const double*, const char*, const int*, const char*, const int*, const double*, double*, const int*, double*, const int*, int*, int*, double*, double*, const int*, int*);
}

static void matvec(int n, double* x, double* y)
{
	y[0]= 2*x[0]-x[1];
	y[1]=-x[0]+2*x[1]-x[2];
	y[2]=-x[1]+2*x[2]-x[3];
	y[3]=-x[2]+2*x[3];
}


int main(void)
{
	//data
	const int n = 4;
	const int nev = 3;
	const int ncv = 4;
	const int lworkl = ncv * ncv + 8 * ncv;

	const double tol = 1.00e-10;
	int ido = 0, info = 0, iparam[11], ipntr[11];
	double resid[n], v[n * n], workd[3 * n], workl[lworkl];

	//setup
	iparam[0] = 1;
	iparam[6] = 1;
	iparam[2] = 300;

	//Arnoldi
	while(true)
	{
		//decomposition
		dsaupd_(&ido, "I", &n, "SA", &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
		//check
		if(ido == 99) break;
		//compute
		if(ido == -1 || ido == 1) matvec(n, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
	}

	if(info != 0)
	{
		std::cout <<"ARPACK failed " << info << "\n";
		return 1;
	}

	//data
	bool select[ncv];
	bool rvec = false;
	double sigma = 0, d[nev], z[n * nev];
	dseupd_(&rvec, "A", select, d, z, &n, &sigma, "I", &n, "SA", &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);

	std::cout <<"Smallest eigenvalue = "<<d[0]<<std::endl;
	std::cout <<"Smallest eigenvalue = "<<d[1]<<std::endl;
	std::cout <<"Smallest eigenvalue = "<<d[2]<<std::endl;
}