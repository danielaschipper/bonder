#pragma once
#include "stdafx.h"
#include "point.h"
#include "readwfn.h"

//here be dragons

//weird mathy stuff 
//this pile of math is used for calculating the intracate math behind wave funtions
class analysisBatch
{
public:
	analysisBatch(wfnData input);
	~analysisBatch();
	//called by the flood fill ajusted from the grid to floating space
	double vauleAtPoint(point Point, void *other);
	//sets up the x,y,z offsets
	void setUpBatch(double x, double y, double z, double res);
	double kinEnergyDensity(double x, double y, double z);
	double hamEnergyDensity(double x, double y, double z);
	double totEnergyDensity(double x, double y, double z);
	double potEnergyDensity(double x, double y, double z);
	double elf(double x, double y, double z);
	//reduced density gradent
	double RDG(double x, double y, double z);
	//arkinsoms energy densities
	double* AKinEng(double x, double y, double z);
	//wavefunction energy densities
	double* WKinEng(double x, double y, double z);
	//reduced energy gradent and the signed vaule of rho
	double RDG_rho(double x, double y, double z, double *srho);

	double Information(double x, double y, double z);
	double LOL(double x, double y, double z);
	double goshEntropy(double x, double y, double z);
	double fisher(double x, double y, double z);

	//x,y,z of atom i
	double atomx(int i);
	double atomy(int i);
	double atomz(int i);
private:
	//clears the vectors
	void vectorReset(double x, double y, double z);
	//calculates the second derivite of each moleculer orbatal at a point
	inline void wavefunctionSecondDerivitivePrimitive(int i);
	void wavefunctionSecondDerivitive();
	//calculates the hessian of each molecular obatal at a point
	inline void wavefuntionHessianPrimitive(int i);
	void wavefuntionhessian();
	//calculates the derivite of each molecular orbatal at a point
	inline void wavefuntionDerivitivePrimitive(int i);
	void wavefuntionDerivitive();
	//caculates the vaule of each moleculer orbatal at a point
	inline void wavefuntionVaulePrimitive(int i);
	void wavefuntionVaule();

	//gets the hessian of the electron density function at a point
	void electronDensityHessian();
	//magic
	double getSignOfSecondEiganVaule();


	double *__restrict__ moWavefuntionDX, *__restrict__ moWavefuntionDy, *__restrict__ moWavefuntionDZ, *__restrict__ distanceFromCenter, *__restrict__ dx, *__restrict__ dy, *__restrict__ dz;
	double *__restrict__ moWavefuntionDXX, *__restrict__ moWavefuntionDYY, *__restrict__ moWavefuntionDZZ, *__restrict__ moWavefuntionVaule, *__restrict__ moWavefuntionHessian, *__restrict__ elecHess;
	//values for basis functions
	double *__restrict__ basisX,*__restrict__ basisY,*__restrict__ basisZ;


	double *__restrict__ molecularOcupancyNumber, *__restrict__ *__restrict__ moleculerOrbatalCoefecents, *__restrict__ centerXvaule, *__restrict__ centerYvaule, *__restrict__ centerZvaule, *__restrict__ primitiveExponantationVaule;
	int *__restrict__ primitiveCenter, *__restrict__ primitiveType;
	double offsetx, offsety, offsetz, res;
	int prims, nmo, centers;
};



