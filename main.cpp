#include "stdafx.h"
#include "analize.h"
#include "readwfn.h"
#include "iface.h"
#include "output.h"
#include <string>
#include <string.h>
#include <iostream>
#include <thread>
#include <boost/asio/io_service.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>


wfnData* init(std::string file)
{
	wfnData *inputFile = readFile(file);
	return inputFile;
}

void drawline(int a, int b, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,bool makeCube)
{
	analysisBatch* batch = new analysisBatch(*inputFile);
	double lowX = (*batch).atomx(a);
	double lowY = (*batch).atomy(a);
	double lowZ = (*batch).atomz(a);

	double highX = (*batch).atomx(b);
	double highY = (*batch).atomy(b);
	double highZ = (*batch).atomz(b);
	
	printf("testing line between %d %d\n",a,b);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 80)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;

	int reps = 1 / jumpScaler;
	bool sucsess = false;
	int flips = 1;
	for (size_t i = 0; i < reps; i++)
	{
		int k = reps / 2 + flips * i / 2;
		flips *= -1;
		double mesured = (*batch).RDG(lowX + k*dx, lowY + k*dy, lowZ + k*dz);
		if (mesured <= cutoff)
		{
			
			analysis analize = analysis();
			analize.setUpAnalysisBatch(lowX + k*dx, lowY + k*dy, lowZ + k*dz, res,batch);

			analize.anilizePoint(0, 0, 0, 0, size, size, cutoff, &sucsess, inputFile, outputfile, batch,makeCube);
			break;
		}
	}
	

	delete batch;
}

struct pdrawArgs
{
	int a; int b; double res; double cutoff; std::string outputfile; int size; wfnData* inputFile; bool makeCube;
	pdrawArgs(int A, int B, double Res, double Cutoff, std::string Outputfile, int Size, wfnData* InputFile,bool MakeCube)
	{
		makeCube = MakeCube;
		a = A;
		b = B;
		res = Res;
		cutoff = Cutoff;
		outputfile = Outputfile;
		size = Size;
		inputFile = InputFile;
	}
};

void pDrawline(void *input)
{
	pdrawArgs* data = (pdrawArgs*)input;
	drawline((*data).a, (*data).b, (*data).res, (*data).cutoff, (*data).outputfile, (*data).size, (*data).inputFile, (*data).makeCube);
	//pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
	if (argc == 1)
	{
		printf("bonder h for help\n");
		return 0;
	}
		

	if (argv[1][0] == 'h')
	{
		printf("the first letter determins the wht the program will do \n p looks at a point and determins the volume around it\n l looks at a line between two atoms\n a looks for all interactions\n g prints out a grid\n h displays this message\n for more detail on an operation type bonder letter");
		return 0;
	}
	
	wfnData *inputFile = 0;
	if (argc != 2)
		inputFile = init(argv[2]);
	std::cout << "data read" << std::endl;
	const int size = 600;
	//letter file x y z res cutoff
	if (argv[1][0] == 'p')
	{
		if (argc != 9 && argc != 10 )
		{
			printf("arguments are bonder p inputFile x y z res cutoff outputFile\n");
			return 0;
		}
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		analize.setUpAnalysisBatch( std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]),batch);
		printf("%f \n", (*batch).RDG(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5])));
		bool sucsess;
		if (argc == 10)
		{
			analize.anilizePoint(0, 0, 0, 0, size, size, std::stod(argv[7]), &sucsess, inputFile, argv[8], batch, !strcmp(argv[9], "true"));
		}
		else
		{
			analize.anilizePoint(0, 0, 0, 0, size, size, std::stod(argv[7]), &sucsess, inputFile, argv[8], batch, true);
		}
		
		if (sucsess)
		{
			printf("point given is in region\n");
		}
		else
		{
			printf("point given is not in region\n");
		}
		return 0;
	}

	//letter file 1 2 res cutoff
	if (argv[1][0] == 'l')
	{
		if (!(argc == 8 || argc == 9))
		{
			printf("arguments are bonder l inputFile atom1 atom2 res cutoff outputFile\n");
			return 0;
		}
		analysisBatch* batch = new analysisBatch(*inputFile);
		if (argc == 8)
			drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], size, inputFile,true);
		else
			drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], size, inputFile, !strcmp(argv[8], "true"));
		return 0;

	}

	//letter file res cutoff output
	if (argv[1][0] == 'a')
	{
		if (argc != 6 && argc != 7)
		{
			printf("arguments are bonder a inputFile res cutoff outputFile\n");
			return 0;
		}
		

		//set up threadpool
		boost::asio::io_service ioService;
		boost::thread_group threadpool;
		boost::asio::io_service::work work(ioService);
		int numOfThreads = std::thread::hardware_concurrency();
		std::cout << numOfThreads << std::endl;
		//numOfThreads = numOfThreads ? numOfThreads : 4;
		numOfThreads = 4;
		
		for(int i = 0; i< numOfThreads; i++)
			threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));


		for (size_t i = 0; i < (*inputFile).nuc; i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				pdrawArgs *lineData;
				if (argc == 7)
				{
					//std::cout << argv[6] << std::endl;
					lineData = new pdrawArgs(i, j, std::stod(argv[3]), std::stod(argv[4]), std::string(argv[5]), size, inputFile, !strcmp(argv[6], "true"));
				}
				else
				{
					lineData = new pdrawArgs(i, j, std::stod(argv[3]), std::stod(argv[4]), std::string(argv[5]) , size, inputFile, true);
				}	
				ioService.post(boost::bind(pDrawline, (void *)lineData));

			}
		}
		std::this_thread::sleep_for(std::chrono::seconds(120));	
		ioService.stop();
		threadpool.join_all();

		return 0;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (argv[1][0] == 'g')
	{
		if (argc != 11 && argc != 12)
		{
			printf("arguments are bonder g inputFile minx miny minz maxx maxy maxz res outputFile\n");
			return 0;
		}
		analysisBatch* batch = new analysisBatch(*inputFile);
		if (argc == 12)
		{
			outputCube(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), argv[10], *inputFile, 1.0, batch, !strcmp(argv[11], "true"));
		}
		else
		{

			outputCube(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), argv[10], *inputFile, 1.0, batch, true);
		}
		
		printf("done");
		return 0;
	}

	//bool sucsess;
	//anilizePoint(0, 0, 0, 0, 2000, 2000, 2501, &sucsess);
	//printf("%d", sucsess);
	printf("bonder h for help\n");
	return 0;
}

/*
int main()
{
	char *args[8] = { "","a", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\h2o-h2o.wfn","0.02","0.5","C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\output" };
	//char *args[9] = { "", "l", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\ethandiol.wfn","7","10", "0.02", "0.5", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\output" };
	//char *args[9] = { "", "p", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\input.wfn","-2.87","1.83","-0.7", "0.02", "0.5", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\output" };
	//char *args[11] = { "", "g", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\ethandiol.wfn","-0.5","-0.5","-0.5","0.5","0.5","0.5", "0.02", "C:\\Users\\ds143\\Documents\\Visual Studio 2015\\Projects\\hydrogen bond project\\Debug\\output" };
	main2(6, args);
	char quit;
	scanf_s("%c", &quit);
	return 0;
}
*/
