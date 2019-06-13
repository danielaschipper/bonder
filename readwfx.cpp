#include "readwfn.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>


using namespace std;

void xspilt(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (!item.empty())
			elems.push_back(item);
	}
}


std::vector<std::string> xspilt(const std::string &s, char delim) {
	std::vector<std::string> elems;
	xspilt(s, delim, elems);
	return elems;
}


void nucliiCount( wfnData* output,ifstream* inputFile)
{
	string line;
	std::getline(*inputFile, line);
	output->nuc = stoi(line);
	std::getline(*inputFile, line);
	output->x = new double[output->nuc];
	output->y = new double[output->nuc];
	output->z = new double[output->nuc];
	output->name = new string[output->nuc];
	output->charge = new int[output->nuc];
}
void MOcount(wfnData* output,ifstream* inputFile)
{
	string line;
	std::getline(*inputFile, line);
	output->MO = stoi(line);
	std::getline(*inputFile, line);
	output->molecularOrbitalOcupancyNumber = new double[output->MO];
}
void atomicNumbers(wfnData* output,ifstream* inputFile)
{
	int i=0;
	string line;
	while (i != output->nuc)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			output->charge[i] = stoi(token);
			i++;
		}
	}
	std::getline(*inputFile, line);
}

void AtomXYZ(wfnData* output,ifstream* inputFile)
{
	int i=0;
	string line;
	double* data = new double[output->nuc * 3];
	while (i != (output->nuc) * 3)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			data[i] = stod(token);
			i++;
		}
	}
	std::getline(*inputFile, line);
	for(i =0; i < output->nuc; i++)
	{
		output->x[i] = data[3*i];
		output->y[i] = data[3*i+1];
		output->z[i] = data[3*i+2];
	}
	delete data;

}

void primCount(wfnData* output,ifstream* inputFile)
{
	string line;
	std::getline(*inputFile, line);
	output->prim = stoi(line);
	std::getline(*inputFile, line);
	output->primitiveCenters = new int[output->prim];
	output->primitiveOrbatalTypes = new int[output->prim];
	output->primitiveExponents = new double[output->prim];
	cout << output->prim << endl;
}

void primCent(wfnData* output,ifstream* inputFile)
{
	int i=0;
	string line;
	while (i != output->prim)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			output->primitiveCenters[i] = stoi(token);
			i++;
		}
	}
	std::getline(*inputFile, line);

}

void primType(wfnData* output,ifstream* inputFile)
{
	int i=0;
	string line;
	while (i != output->prim)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			output->primitiveOrbatalTypes[i] = stoi(token);
			i++;
		}
	}
	std::getline(*inputFile, line);

}

void primExp(wfnData* output,ifstream* inputFile)
{
	int i=0;
	string line;
	while (i != output->prim)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			output->primitiveExponents[i] = stod(token);
			i++;
		}
	}
	std::getline(*inputFile, line);

}

void moo(wfnData* output,ifstream* inputFile)
{
	int i =0;
	string line;
	while (i != output->MO)
	{
		std::getline(*inputFile, line);
		vector<string> tokens = xspilt(line,' ');
		for (string token: tokens)
		{
			output->molecularOrbitalOcupancyNumber[i] = stod(token);
			i++;
		}
	}
	std::getline(*inputFile, line);

}

void mos(wfnData* output,ifstream* inputFile)
{
	string line;
	(*output).molecularOrbatalCoeficents = new double*[(*output).MO];
	for (int i = 0; i < (*output).prim; i++)
        {
                (*output).molecularOrbatalCoeficents[i] = new double[(*output).prim];
        }

	for (int j = 0; j < output->MO; j++ )
	{
		std::getline(*inputFile, line);
		std::getline(*inputFile, line);
		std::getline(*inputFile, line);
		int i = 0;
		while (i != output->prim)
		{
			std::getline(*inputFile, line);
			vector<string> tokens = xspilt(line,' ');
			for (string token: tokens)
			{
				output->molecularOrbatalCoeficents[j][i] = stod(token);
				i++;
			}
		}

	}
	std::getline(*inputFile, line);
}

wfnData* readWfx(string file)
{
	wfnData* output = new wfnData;
	ifstream inputFile(file);
	if (!inputFile)
	{
		printf("wavefunction file not found\n");
		exit(1);
	}


	string line;
	while(std::getline(inputFile, line))
	{
		if(line[0] == '<')
		{
			if(line == "<Number of Nuclei>")
				nucliiCount(output,&inputFile);
			else if (line == "<Number of Occupied Molecular Orbitals>")
				MOcount(output,&inputFile);
			else if (line == "<Atomic Numbers>")
				atomicNumbers(output,&inputFile);
			else if (line == "<Nuclear Cartesian Coordinates>")
				AtomXYZ(output,&inputFile);
			else if (line == "<Number of Primitives>")
				primCount(output,&inputFile);
			else if (line == "<Primitive Centers>")
				primCent(output,&inputFile);
			else if (line == "<Primitive Types>")
				primType(output,&inputFile);
			else if (line == "<Primitive Exponents>")
				primExp(output,&inputFile);
			else if (line == "<Molecular Orbital Occupation Numbers>")
				moo(output,&inputFile);
			else if (line == "<Molecular Orbital Primitive Coefficients>")
				mos(output,&inputFile);




		}
	}
	return output;

}
