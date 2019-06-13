#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "readwfn.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <math.h>

using namespace std;

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (!item.empty())
			elems.push_back(item);
	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

//decode fortran doubles
double DFD(string input)
{
	std::vector<std::string> seglist = split(input,'D');
	double mines = stod(seglist[0]);
	double power = stoi(seglist[1]);
	return mines * pow(10, power);

}


wfnData* readFile(string file)
{
	
	wfnData* output = new wfnData;
	ifstream inputFile(file);
	
	if (!inputFile)
	{
		printf("wavefunction file not found\n");
		exit(1);
	}

	string line;
	std::getline(inputFile, line);
	
	//read first line
	//token 2 is no mo
	//token 4 is no prim
	//token 6 is no nuci
	std::getline(inputFile, line);
	vector<string> tokens = split(line,' ');
	(*output).MO = stoi(tokens[1]);
	(*output).prim = stoi(tokens[4]);
	(*output).nuc = stoi(tokens[6]);
	tokens.clear();
	//set up for nucli
	(*output).x = new double[(*output).nuc];
	(*output).y = new double[(*output).nuc];
	(*output).z = new double[(*output).nuc];
	(*output).charge = new int[(*output).nuc];
	(*output).name = new string[(*output).nuc];

	//set up prims
	(*output).primitiveCenters = new int[(*output).prim];
	(*output).primitiveOrbatalTypes = new int[(*output).prim];
	(*output).primitiveExponents = new double[(*output).prim];

	//set up MO
	(*output).molecularOrbitalOcupancyNumber = new double[(*output).MO];
	(*output).MOType = new double[(*output).MO];
	(*output).MOeng = new double[(*output).MO];
	(*output).molecularOrbatalCoeficents = new double*[(*output).MO];
	for (int i = 0; i < (*output).prim; i++)
	{
		(*output).molecularOrbatalCoeficents[i] = new double[(*output).prim];
	}

	//read atomic data
	for (int i = 0; i < (*output).nuc; i++)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		//token 1 is letter, 4 is x, 5 is y, 6 is z, 9 is charge
		(*output).name[i] = tokens[0];
		(*output).x[i] = stod(tokens[4]);
		(*output).y[i] = stod(tokens[5]);
		(*output).z[i] = stod(tokens[6]);
		(*output).charge[i] = stoi(tokens[9]);
		tokens.clear();
	}
	int noOfLines = (*output).prim / 20, LeftOver = (*output).prim % 20;
	//read center assinments
	for (int i = 0; i < noOfLines; i++)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 2; j < 22; j++)
		{
			if (tokens[j].length() > 3){
				string mess = tokens[j];
				if (mess.length() % 3 == 2)
				{
					int num = stoi(mess.substr(0,2));
					(*output).primitiveCenters[i * 20 + j - 2] = num;
					mess = mess.erase(0, 2);
					j++;
				}
				for (; j < 22; j++)
				{
					(*output).primitiveCenters[i * 20 + j - 2] = stoi(mess.substr(0,3));
					mess = mess.erase(0, 3);
				}
				
			}
			else
			{
				(*output).primitiveCenters[i * 20 + j -2] = stoi(tokens[j]);
				//printf("%d %d   ", (*output).primitiveCenters[i * 20 + j - 2],(i * 20 + j - 2));
			}
		}
		tokens.clear();
		//printf("\n");
	}
	//read left over center assinments
	if (LeftOver != 0)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 2; j < LeftOver + 2; j++)
		{
			if (tokens[j].length() > 3){
				string mess = tokens[j];
				if (mess.length() % 3 == 2)
				{
					(*output).primitiveCenters[noOfLines * 20 + j - 2] = stoi(mess.substr(0,2));
					mess = mess.erase(0, 2);
					j++;
				}
				for (; j < LeftOver +2; j++)
				{
					(*output).primitiveCenters[noOfLines * 20 + j - 2] = stoi(mess.substr(0,3));
					mess = mess.erase(0, 3);
					j++;
				}

			}
			else
				(*output).primitiveCenters[noOfLines * 20 + j - 2] = stoi(tokens[j]);
			//printf("%d ", (*output).centerAssinments[noOfLines * 20 + j - 2]);
		}
		tokens.clear();
		//printf("\n");
	}

	//read type assinments
	for (int i = 0; i < noOfLines; i++)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 2; j < 22; j++)
		{
			(*output).primitiveOrbatalTypes[i * 20 + j - 2] = stoi(tokens[j]);
			//printf("%d ", (*output).typeAssinments[i * 20 + j - 2]);
		}
		tokens.clear();
	}

	//read left over type assinments
	if (LeftOver != 0)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 2; j < LeftOver + 2; j++)
		{
			(*output).primitiveOrbatalTypes[noOfLines * 20 + j - 2] = stoi(tokens[j]);
			//printf("%d ", (*output).typeAssinments[noOfLines * 20 + j - 2]);
		}
		tokens.clear();
		//printf("\n");
	}

	noOfLines = (*output).prim / 5;
	LeftOver = (*output).prim % 5;

	//read exponents
	for (int i = 0; i < noOfLines; i++)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 1; j < 6; j++)
		{
			(*output).primitiveExponents[i * 5 + j - 1] = DFD(tokens[j]);
			//printf("%f ", (*output).exponents[i * 5 + j - 1]);
		}
		tokens.clear();
	}

	//read left over exponents
	if (LeftOver != 0)
	{
		std::getline(inputFile, line);
		tokens = split(line, ' ');
		for (int j = 1; j < LeftOver + 1; j++)
		{
			(*output).primitiveExponents[noOfLines * 5 + j - 1] = DFD(tokens[j]);
			//printf("%f ", (*output).exponents[noOfLines * 5 + j - 1]);
		}
		tokens.clear();
	}
	//read MO
	for (int k = 0; k < (*output).MO; k++)
	{
		//read mo inintalisation
		//token 4 is type, 8 is occnum, 12 is energy
		std::getline(inputFile, line);
		while (line[0] == 'C'||line[0] == 'T'||line[0] == 'E')
			std::getline(inputFile, line);
		tokens = split(line, ' ');
		//cout << line;
		(*output).MOType[k] = stod(tokens[3]);
		(*output).molecularOrbitalOcupancyNumber[k] = stod(tokens[7]);
		(*output).MOeng[k] = stod(tokens[11]);
		//printf("%f \n", (*output).MOOOC[k]);
		tokens.clear();
		// read mo coefecents
		for (int i = 0; i < noOfLines; i++)
		{
			std::getline(inputFile, line);
			tokens = split(line, ' ');
			for (int j = 0; j < 5; j++)
			{
				(*output).molecularOrbatalCoeficents[k][i * 5 + j] = DFD(tokens[j]);
				//printf("%f ", (*output).MOCO[j][i * 5 + j]);
			}
			tokens.clear();
			//printf("\n\n");
		}

		//read left over MO coefecents
		if (LeftOver != 0)
		{
			std::getline(inputFile, line);
			tokens = split(line, ' ');
			for (int j = 0; j < LeftOver; j++)
			{
				(*output).molecularOrbatalCoeficents[k][noOfLines * 5 + j] = DFD(tokens[j]);
				//printf("%f ", (*output).MOCO[j][noOfLines * 5 + j]);
			}
			tokens.clear();
			//printf("\n\n");
		}
		//printf("\n\n");
	}
	
	inputFile.close();
	return output;
}
