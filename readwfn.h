#pragma once
#include <string>
#include <iostream>

//this class reads in .wfn files
struct wfnData
{
	int prim, nuc, MO;
	
	//atoms
	double *x, *y, *z ;
	std::string *name;
	int *charge;

	//primitives
	int* primitiveCenters, *primitiveOrbatalTypes;
	double* primitiveExponents;

	//MO
	double *MOType, *molecularOrbitalOcupancyNumber;
	double *MOeng, **molecularOrbatalCoeficents;

	void print()
	{
		for(int h = 0; h < MO; h++)
		{
			std::cout<< molecularOrbitalOcupancyNumber[h] << std::endl;
		}

	}
};

wfnData* readFile(std::string file);
wfnData* readWfx(std::string file);
