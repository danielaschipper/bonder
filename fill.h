#pragma once

#include "stdafx.h"
#include "LLI.h"
#include "iface.h"
#include "point.h"


struct gridPoint
{
	
	LLE* edges;
	LLi* internalPoints;
	
	
	gridPoint()
	{
		edges = new LLE;
		internalPoints = new LLi;
	}
};

struct grid
{
	gridPoint *data;
	int x, y;
};

//runs the flood fil algorytem
grid fill(int x, int y, int z, void * other, int Xsize, int Ysize, double cutOff, bool *sucsess, analysisBatch* batch);
gridPoint* getPoint(grid* Grid, int x, int y);
