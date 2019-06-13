#include "stdafx.h"
#include "iface.h"
#include "fill.h"
#include <math.h>

#define PI           3.14159265358979323846
#define abso(x) ((x > 0)? x: -x)
#define A(x,y) (elecHess[y * 3 +x - 4])

const double rhoCutoff = 0.1;

//must be set up


int type2ix[56] = { 0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5 };
int type2iy[56] = { 0, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 };
int type2iz[56] = { 0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 1, 0, 2, 2, 1, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0, 5, 4, 3, 2, 1, 0, 4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0 };


//exponental function
inline double power(double b, int p)
{
	double result = 1;
	for (; p; --p)
	{
		result *= b;
	}
	return result;
}

//setuo
void analysisBatch::setUpBatch(double x, double y, double z, double resalution)
{
	offsetx = x;
	offsety = y;
	offsetz = z;
	res = resalution;
}

//more setup
analysisBatch::analysisBatch(wfnData input)
{
	nmo = input.MO;
	moWavefuntionDX = new double[nmo];
	moWavefuntionDy = new double[nmo];
	moWavefuntionDZ = new double[nmo];
	moWavefuntionDXX = new double[nmo];
	moWavefuntionDYY = new double[nmo];
	moWavefuntionDZZ = new double[nmo];
	moWavefuntionVaule = new double[nmo];
	moWavefuntionHessian = new double[nmo * 9];
	
	elecHess = new double[9];

	centers = input.nuc;
	
	prims = input.prim;
	basisX = new double[prims];
	basisY = new double[prims];
	basisZ = new double[prims];

	distanceFromCenter = new double[centers];
	dx = new double[centers];
	dy = new double[centers];
	dz = new double[centers];


	molecularOcupancyNumber = input.molecularOrbitalOcupancyNumber;
	moleculerOrbatalCoefecents = input.molecularOrbatalCoeficents;
	

	primitiveCenter = input.primitiveCenters;
	primitiveType = input.primitiveOrbatalTypes;

	primitiveExponantationVaule = input.primitiveExponents;

	centerXvaule = input.x;
	centerYvaule = input.y;
	centerZvaule = input.z;




}

//cleanup
analysisBatch::~analysisBatch()
{
	delete moWavefuntionDX;
	delete moWavefuntionDy;
	delete moWavefuntionDZ;
	delete moWavefuntionDXX;
	delete moWavefuntionDYY;
	delete moWavefuntionDZZ;
	delete moWavefuntionVaule;
	delete moWavefuntionHessian;
	delete elecHess;
	delete distanceFromCenter;
	delete dx;
	delete dy;
	delete dz;
}


//clears MO wavefunction data
void analysisBatch::vectorReset(double x,double y, double z)
{
	z += 0;
	for (int i = 0; i < centers; ++i)
	{
		dx[i] = (x - centerXvaule[i]);
		dy[i] = (y - centerYvaule[i]);
		dz[i] = (z - centerZvaule[i]);

		distanceFromCenter[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	}
	for (int i = 0; i < nmo; ++i)
	{
		moWavefuntionDX[i] = 0;
		moWavefuntionDy[i] = 0;
		moWavefuntionDZ[i] = 0;
		moWavefuntionDXX[i] = 0;
		moWavefuntionDYY[i] = 0;
		moWavefuntionDZZ[i] = 0;
		moWavefuntionVaule[i] = 0;
		#pragma unroll
		for (int j = 0; j < 9; ++j)
		{
			moWavefuntionHessian[i * 9 + j] = 0;
		}
	}
	z += 0;
}

//calculates MO wave function data fir each primitive
inline void analysisBatch::wavefunctionSecondDerivitivePrimitive(int i)
{
	int ix, iy, iz,center;
	double tx, ty, tz, GTFdx, GTFdy, GTFdz;
	double expt;
	center = primitiveCenter[i] - 1;

	if (distanceFromCenter[center] * primitiveExponantationVaule[i] > 15)
		return;

	ix = type2ix[primitiveType[i] - 1];
	iy = type2iy[primitiveType[i] - 1];
	iz = type2iz[primitiveType[i] - 1];


	expt = exp(-distanceFromCenter[center] * primitiveExponantationVaule[i]);
	tx = (ix > 1) ? ix * (ix - 1) * power(dx[center], (ix - 2)) : 0;
	ty = (iy > 1) ? iy * (iy - 1) * power(dy[center], (iy - 2)) : 0;
	tz = (iz > 1) ? iz * (iz - 1) * power(dz[center], (iz - 2)) : 0;

	basisX[i] = power(dy[center], iy) * power(dz[center], iz) *expt * (tx + 2 * primitiveExponantationVaule[i] * power(dx[center], ix)*(-2 * ix + 2 * primitiveExponantationVaule[i] * dx[center] * dx[center] - 1));
	basisY[i] = power(dz[center], iz) * power(dx[center], ix) *expt * (ty + 2 * primitiveExponantationVaule[i] * power(dy[center], iy)*(-2 * iy + 2 * primitiveExponantationVaule[i] * dy[center] * dy[center] - 1));
	basisZ[i] = power(dx[center], ix) * power(dy[center], iy) *expt * (tz + 2 * primitiveExponantationVaule[i] * power(dz[center], iz)*(-2 * iz + 2 * primitiveExponantationVaule[i] * dz[center] * dz[center] - 1));

}

void analysisBatch::wavefunctionSecondDerivitive()
{
	//getWaveFnVaule(x, y, z);
	for (int i = 0; i < prims; ++i)
	{
		basisX[i] = 0;
		basisY[i] = 0;
		basisZ[i] = 0;
	}
	#pragma omp parallel for 	
	for (int i = 0; i < prims; ++i)
	{
		wavefunctionSecondDerivitivePrimitive(i);
	}
	//TODO matrix mult
	#pragma omp parallel for
	for (int j = 0; j < nmo; ++j)
	{
		for (int i = 0; i < prims; ++i)
		{
			moWavefuntionDXX[j] += moleculerOrbatalCoefecents[j][i] * basisX[i];
			moWavefuntionDYY[j] += moleculerOrbatalCoefecents[j][i] * basisY[i];
			moWavefuntionDZZ[j] += moleculerOrbatalCoefecents[j][i] * basisZ[i];
		}
	}


}


//computes the hessian matrix for the MO wavefunction primitives
inline void analysisBatch::wavefuntionHessianPrimitive(int i)
{
	int ix, iy, iz;
	double GTFdxy, GTFdyz, GTFdzx,ttx,tty,ttz;
	double expt;
	int center = primitiveCenter[i] - 1;

	if (distanceFromCenter[center] * primitiveExponantationVaule[i] > 15)
		return;

	ix = type2ix[primitiveType[i] - 1];
	iy = type2iy[primitiveType[i] - 1];
	iz = type2iz[primitiveType[i] - 1];


	expt = exp(-distanceFromCenter[center] * primitiveExponantationVaule[i]);


	ttx = ((ix > 1) ? ix * (ix - 1) * power(dx[center], (ix - 2)) : 0) - 2 * primitiveExponantationVaule[i] * power(dx[center], ix + 1);
	tty = ((iy > 1) ? iy * (iy - 1) * power(dy[center], (iy - 2)) : 0) - 2 * primitiveExponantationVaule[i] * power(dy[center], iy + 1);
	ttz = ((iz > 1) ? iz * (iz - 1) * power(dz[center], (iz - 2)) : 0) - 2 * primitiveExponantationVaule[i] * power(dz[center], iz + 1);



	basisX[i] = GTFdxy = power(dz[center], iz) * expt * ttx * tty;
	basisY[i] = GTFdyz = power(dx[center], ix) * expt * ttz * tty;
	basisZ[i] = GTFdzx = power(dy[center], iy) * expt * ttx * ttz;


}

void analysisBatch::wavefuntionhessian()
{
	wavefunctionSecondDerivitive();
	for (int i = 0; i < prims; ++i)
	{
		basisX[i] = 0;
		basisY[i] = 0;
		basisZ[i] = 0;
	}
	//getWaveFnVaule(x, y, z);	
	#pragma omp parallel for 
	for (int i = 0; i < prims; ++i)
	{
		wavefuntionHessianPrimitive(i);
	}
	//TODO matrix

	#pragma omp parallel for
	for (int j = 0; j < nmo; ++j)
	{
		for (int i = 0; i < prims; ++i)
		{
			moWavefuntionHessian[j * 9 + 1] += moleculerOrbatalCoefecents[j][i] * basisX[i];
			moWavefuntionHessian[j * 9 + 2] += moleculerOrbatalCoefecents[j][i] * basisY[i];
			moWavefuntionHessian[j * 9 + 5] += moleculerOrbatalCoefecents[j][i] * basisZ[i];
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < nmo; ++i)
	{
		moWavefuntionHessian[i * 9 ] = moWavefuntionDXX[i];
		//moWavefuntionHessian[i * 9 + 3] = moWavefuntionHessian[i * 9 + 1];
		moWavefuntionHessian[i * 9 + 4] = moWavefuntionDYY[i];
		//moWavefuntionHessian[i * 9 + 6] = moWavefuntionHessian[i * 9 + 2];
		//moWavefuntionHessian[i * 9 + 7] = moWavefuntionHessian[i * 9 + 5];
		moWavefuntionHessian[i * 9 + 8] = moWavefuntionDZZ[i];
	}
}

//calculates the firts derivitive of the MO wavefunction
inline void analysisBatch::wavefuntionDerivitivePrimitive(int i)
{
	int ix, iy, iz;
	double tx, ty, tz, GTFdx, GTFdy, GTFdz;
	double expt;
	int center = primitiveCenter[i] - 1;

	if (distanceFromCenter[center] * primitiveExponantationVaule[i] > 15)
		return;

	ix = type2ix[primitiveType[i] - 1];
	iy = type2iy[primitiveType[i] - 1];
	iz = type2iz[primitiveType[i] - 1];

	expt = exp(-distanceFromCenter[center] * primitiveExponantationVaule[i]);
	tx = ix ? ix * power(dx[center], (ix - 1)) : 0;
	ty = iy ? iy * power(dy[center], (iy - 1)) : 0;
	tz = iz ? iz * power(dz[center], (iz - 1)) : 0;

	basisX[i] = power(dy[center], iy) * power(dz[center], iz) * expt * (tx - 2 * primitiveExponantationVaule[i] * power(dx[center], ix + 1));
	basisY[i] = power(dx[center], ix) * power(dz[center], iz) * expt * (ty - 2 * primitiveExponantationVaule[i] * power(dy[center], iy + 1));
	basisZ[i] = power(dx[center], ix) * power(dy[center], iy) * expt * (tz - 2 * primitiveExponantationVaule[i] * power(dz[center], iz + 1));

}

void analysisBatch::wavefuntionDerivitive()
{
	
	for (int i = 0; i < prims; ++i)
	{
		basisX[i] = 0;
		basisY[i] = 0;
		basisZ[i] = 0;
	}
	//getWaveFnVaule(x, y, z);	
	#pragma omp parallel for 
	for (int i = 0; i < prims; ++i)
	{
		wavefuntionDerivitivePrimitive(i);
	}
	#pragma omp parallel for
	for (int j = 0; j < nmo; ++j)
	{
		for (int i = 0; i < prims; ++i)
		{
			moWavefuntionDX[j] += moleculerOrbatalCoefecents[j][i] * basisX[i];
			moWavefuntionDy[j] += moleculerOrbatalCoefecents[j][i] * basisY[i];
			moWavefuntionDZ[j] += moleculerOrbatalCoefecents[j][i] * basisZ[i];
		}
	}
}



inline void analysisBatch::wavefuntionVaulePrimitive(int i)
{
	int ix, iy, iz;
	double GTFval;
	double expt;
	int center = primitiveCenter[i] - 1;
	if (distanceFromCenter[center] * primitiveExponantationVaule[i] > 15)
		return;

	int type = primitiveType[i] - 1;
	ix = type2ix[type];
	iy = type2iy[type];
	iz = type2iz[type];

	expt = exp(-distanceFromCenter[center] * primitiveExponantationVaule[i]);
	basisY[i] = expt * power(dx[center], ix) * power(dy[center], iy) * power(dz[center], iz);

}

void analysisBatch::wavefuntionVaule()
{
	//getWaveFnVaule(x, y, z);
	for (int i = 0; i < prims; ++i)
		basisY[i] = 0;	
	#pragma omp parallel for 
	for (int i = 0; i < prims; ++i)
	{
		wavefuntionVaulePrimitive(i);
	}

	#pragma omp parallel for
	for (int j = 0; j < nmo; ++j)
	{
		for (int i = 0; i < prims; ++i)
		{
			moWavefuntionVaule[j] += moleculerOrbatalCoefecents[j][i] * basisY[i];
		}
	}

}

double analysisBatch::atomx(int i)
{
	return centerXvaule[i];
}

double analysisBatch::atomy(int i)
{
	return centerYvaule[i];
}

double analysisBatch::atomz(int i)
{
	return centerZvaule[i];
}




const int devive = 1;
double analysisBatch::kinEnergyDensity(double x, double y, double z)
{
	vectorReset(x, y, z);
	wavefuntionDerivitive();
	double kineng = 0;
	for (int i = 0; i < nmo; ++i)
	{
		kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i]* moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
	}
	return kineng / 2;
	
}

double analysisBatch::hamEnergyDensity(double x, double y, double z)
{
	vectorReset(x, y, z);
	wavefuntionVaule();
	wavefunctionSecondDerivitive();
	double hameng = 0;
	for (int i = 0; i < nmo; ++i)
	{
		hameng += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * (moWavefuntionDXX[i] + moWavefuntionDYY[i] + moWavefuntionDZZ[i]);
	}
	return - hameng / 2;
}

double analysisBatch::totEnergyDensity(double x, double y, double z)
{
	return -hamEnergyDensity(x, y, z);
}

double analysisBatch::potEnergyDensity(double x, double y, double z)
{
	return -(hamEnergyDensity(x, y, z) + kinEnergyDensity(x,y,z));
}

void analysisBatch::electronDensityHessian()
{
#pragma unroll
	for (int i = 0; i < 9; ++i)
	{
		elecHess[i] = 0;
	}

	for (int i = 0; i < nmo; ++i)
	{
		elecHess[0] += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 0]);
		elecHess[1] += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDy[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 1]);
		elecHess[2] += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDZ[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 2]);
		elecHess[4] += molecularOcupancyNumber[i] * (moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 4]);
		elecHess[5] += molecularOcupancyNumber[i] * (moWavefuntionDy[i] * moWavefuntionDZ[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 5]);
		elecHess[8] += molecularOcupancyNumber[i] * (moWavefuntionDZ[i] * moWavefuntionDZ[i] + moWavefuntionVaule[i] * moWavefuntionHessian[9 * i + 8]);
	}
	elecHess[3] = elecHess[1];
	elecHess[6] = elecHess[2];
	elecHess[7] = elecHess[5];

#pragma unroll
	for (int i = 0; i < 9; ++i)
	{
		elecHess[i] *= 2;
	}

}

double analysisBatch::getSignOfSecondEiganVaule()
{
	//get eiganVaules
	double eigan[3],phi;

	double p1 = A(1, 2) * A(1, 2) + A(1, 3) * A(1, 3) + A(2, 3) * A(2, 3);
	if (p1 == 0)
	{
		// A is diagonal.
		eigan[0] = A(1, 1);
		eigan[1] = A(2, 2);
		eigan[2] = A(3, 3);
	}
	else
	{
		double q = (A(1,1) + A(2,2) + A(3,3)) / 3;
		double p2 = (A(1, 1) - q) * (A(1, 1) - q) + (A(2, 2) - q) * (A(2, 2) - q) + (A(3, 3) - q) * (A(3, 3) - q) + 2 * p1;
		double p = sqrt(p2 / 6);
		double B[9];
		for (int i = 0; i < 9; ++i)
		{
			B[i] = (1 / p) * (elecHess[i] - q * ((i % 4) ? 0 : 1)) ; // I is the identity matrix;
		}
		
		double r = (B[0] * B[4] * B[8] - B[0] * B[5] * B[7] - B[1] * B[3] * B[8] + B[1] * B[5] * B[6] + B[2] * B[3] * B[7] - B[2] * B[4] * B[6]) / 2;

		// In exact arithmetic for a symmetric matrix - 1 <= r <= 1;
		// but computation error can leave it slightly outside this range.
		if (r <= -1)
			phi = PI / 3;
		else if (r >= 1)
			phi = 0;
		else
			phi = acos(r) / 3;

		//the eigenvalues satisfy eig3 <= eig2 <= eig1
		eigan[0] = q + 2 * p * cos(phi);
		eigan[2] = q + 2 * p * cos(phi + (2 * PI / 3));
		eigan[1] = 3 * q - eigan[0] - eigan[2];    // since trace(A) = eig1 + eig2 + eig3;
	}
	//printf("%f\n",eigan[1]);
	return ((eigan[1] > 0) ? 1 : -1);
}

double analysisBatch::RDG_rho(double x, double y, double z,double *srho)
{

	vectorReset(x, y, z);
	wavefuntionVaule();
	wavefuntionDerivitive();
	wavefuntionhessian();
	double  Rhox = 0, rhoy = 0, rhoz = 0, rhosq;
	double rho = 0;
	//#pragma omp parallel for reduction(+:rho,Rhox,rhoy,rhoz)
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
		rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
		rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
	}

	Rhox *= 2;
	rhoy *= 2;
	rhoz *= 2;

	rhosq = Rhox*Rhox + rhoy*rhoy + rhoz*rhoz;

	double absrho = (rho < 0) ? -rho : rho;
	electronDensityHessian();
	*srho = absrho * getSignOfSecondEiganVaule();

	if (rho == 0)
	{
		return 999;
	}

	 
	
	if (absrho > rhoCutoff)
		return 999;

	return 0.161620459673995 * sqrt(rhosq) / (pow(rho, 4.0 / 3));
	//0.161620459673995 =  1/(2*(3*pi^2)^(1/3))
}

double analysisBatch::Information(double x, double y, double z)
{
	double rho = 0;
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
	}
	return -rho * log(rho);
}

double analysisBatch::LOL(double x, double y, double z)
{
	double rho = 0;
	double kineng = 0;
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
	}
	kineng /= 2;
	double Fc = 2.871234000; //magic number do not touch
	double Dh = Fc*pow(rho,(5.0 / 3.0));
	kineng = Dh / kineng;
	return 1.0/(1.0 + kineng);
}

double analysisBatch::goshEntropy(double x, double y, double z)
{
	double rho = 0;
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
	}

	double ck = 2.871234; //magic
	double TFkin = ck* pow(rho, (5.0 / 3.0));

	double kineng = 0;
	for (int i = 0; i < nmo; ++i)
	{
		kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
	}

	kineng /= 2;

	double xlap = 0 , ylap = 0, zlap = 0;
	for (int i = 0; i < nmo; ++i)
	{
		xlap += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionVaule[i] * moWavefuntionDXX[i]);
		ylap += molecularOcupancyNumber[i] * (moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionVaule[i] * moWavefuntionDYY[i]);
		ylap += molecularOcupancyNumber[i] * (moWavefuntionDZ[i] * moWavefuntionDZ[i] + moWavefuntionVaule[i] * moWavefuntionDZZ[i]);
	}
	kineng -= (xlap + ylap + zlap) / 4;

	double rlambda = 5.0 / 3.0 + log(4.0*PI*ck / 3.0);
	double rlogterm = log(kineng / TFkin);
	return 1.5*rho*(rlambda + rlogterm);
}

double analysisBatch::fisher(double x, double y, double z)
{
	double rho = 0, Rhox = 0, rhoy = 0, rhoz = 0, rhosq;
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
		rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
		rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
	}

	Rhox *= 2;
	rhoy *= 2;
	rhoz *= 2;

	rhosq = Rhox*Rhox + rhoy*rhoy + rhoz*rhoz;
	return rhosq / rho;
}



double analysisBatch::RDG(double x, double y, double z)
{
	vectorReset(x, y, z);
	wavefuntionVaule();
	wavefuntionDerivitive();
	double rho = 0,Rhox = 0,rhoy = 0, rhoz = 0,rhosq;
	//#pragma omp parallel for reduction(+:rho,Rhox,rhoy,rhoz)
	for (int i = 0; i < nmo; ++i)
	{
		rho  += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
		rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
		rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
	}

	Rhox *= 2;
	rhoy *= 2;
	rhoz *= 2;

	rhosq = Rhox*Rhox + rhoy*rhoy + rhoz*rhoz;

	if (rhosq == 0 || rho == 0)
	{
		return 999;
	}
	double absrho = (rho < 0) ? -rho : rho;
	if (absrho > rhoCutoff)
		return 999;

	return 0.161620459673995 * sqrt(rhosq) / (pow(rho, 4.0/3));
	//0.161620459673995D0 =  1/(2*(3*pi^2)^(1/3))
}

double* analysisBatch::AKinEng(double x, double y, double z)
{
	double * output = new double[3];
	//derivInit(x, y, z);
	//wfnDerv();
	//wfnval();
	//wfnsdv();

	double rho = 0, Rhox = 0, rhoy = 0, rhoz = 0,lapx = 0,lapy = 0, lapz =0;
	//#pragma omp parallel for reduction(+:rho,Rhox,rhoy,rhoz,lapx,lapy,lapz)
	for (int i = 0; i < nmo; ++i)
	{
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
		rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
		rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
		lapx += molecularOcupancyNumber[i] * (power(moWavefuntionDX[i], 2) + moWavefuntionVaule[i] * moWavefuntionDXX[i]);
		lapy += molecularOcupancyNumber[i] * (power(moWavefuntionDy[i], 2) + moWavefuntionVaule[i] * moWavefuntionDYY[i]);
		lapz += molecularOcupancyNumber[i] * (power(moWavefuntionDZ[i], 2) + moWavefuntionVaule[i] * moWavefuntionDZZ[i]);
	}

	Rhox *= 2;
	rhoy *= 2;
	rhoz *= 2;
	lapx *= 2;
	lapy *= 2;
	lapz *= 2;
	double rhosq = Rhox * Rhox + rhoy * rhoy + rhoz * rhoz;
	double lap = lapx + lapy + lapz;
	output[0] = (3.0 * pow(3.0 * PI * PI, 2.0 / 3.0) * pow(rho, 5.0 / 3.0) / 10.0) + rhosq / (72 * rho) + lap / 6.0;
	output[1] = lap / 4.0 - 2 * output[0];
	output[2] = output[0] + output[1];
	return output;
}

double* analysisBatch::WKinEng(double x, double y, double z)
{
	double * output = new double[3];
	output[0] = kinEnergyDensity(x,y,z); 
	output[2] = -hamEnergyDensity(x,y,z);
	output[1] = output[2]-output[0];
}

double analysisBatch::vauleAtPoint(point Point, void * other)
{
	return RDG(Point.x * res + offsetx, Point.y * res + offsety, Point.z * res + offsetz);
}

double analysisBatch::elf(double x, double y, double z)
{
	double rho = 0, Rhox = 0, rhoy = 0, rhoz = 0;
	double kineng = 0;
	double Fc = 2.871234000; //magic number do not touch
	//#pragma omp parallel for reduction(+:rho,Rhox,rhoy,rhoz,kineng)
	for (int i = 0; i < nmo; ++i)
	{
		kineng += molecularOcupancyNumber[i] * (moWavefuntionDX[i] * moWavefuntionDX[i] + moWavefuntionDy[i] * moWavefuntionDy[i] + moWavefuntionDZ[i] * moWavefuntionDZ[i]);
		rho += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionVaule[i];
		Rhox += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDX[i];
		rhoy += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDy[i];
		rhoz += molecularOcupancyNumber[i] * moWavefuntionVaule[i] * moWavefuntionDZ[i];
	}
	Rhox *= 2;
	rhoy *= 2;
	rhoz *= 2;
	double Dh = Fc*pow(rho,(5.0/ 3.0));
	kineng = kineng / 2.0 - (Rhox*Rhox + rhoy*rhoy + rhoz*rhoz) / rho / 8.0;
	return 1 / (1 + (kineng / Dh)*(kineng / Dh));
}
