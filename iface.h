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
	//arkinsoms reduced energy gradent
	double* AKinEng(double x, double y, double z);
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


	double *__restrict moWavefuntionDX, *__restrict moWavefuntionDy, *__restrict moWavefuntionDZ, *__restrict distanceFromCenter, *__restrict dx, *__restrict dy, *__restrict dz;
	double *__restrict moWavefuntionDXX, *__restrict moWavefuntionDYY, *__restrict moWavefuntionDZZ, *__restrict moWavefuntionVaule, *__restrict moWavefuntionHessian, *__restrict elecHess;



	double *__restrict molecularOcupancyNumber, *__restrict *__restrict moleculerOrbatalCoefecents, *__restrict centerXvaule, *__restrict centerYvaule, *__restrict centerZvaule, *__restrict primitiveExponantationVaule;
	int *__restrict primitiveCenter, *__restrict primitiveType;
	double offsetx, offsety, offsetz, res;
	int prims, nmo, centers;
};



