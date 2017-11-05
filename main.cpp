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
#include <fstream>


const int SIZE = 600;

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

	printf("testing line between %d and %d\n",a,b);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	printf("%d\n",reps);
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

void drawtrig(int a, int b,int c, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,bool makeCube)
{

	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = (*batch).atomx(c);
	double lowY = (*batch).atomy(c);
	double lowZ = (*batch).atomz(c);

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	printf("%d\n",reps);
	for (size_t i = 0; i < reps; i++)
	{
		int k = i;
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

void drawquad(int a, int b,int c,int d, double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,bool makeCube)
{

	analysisBatch* batch = new analysisBatch(*inputFile);
	double highX = ((*batch).atomx(a) + (*batch).atomx(b))/2;
	double highY = ((*batch).atomy(a) + (*batch).atomy(b))/2;
	double highZ = ((*batch).atomz(a) + (*batch).atomz(b))/2;

	double lowX = ((*batch).atomx(c)+(*batch).atomx(d))/2;
	double lowY = ((*batch).atomy(c)+(*batch).atomx(d))/2;
	double lowZ = ((*batch).atomz(c)+(*batch).atomx(d))/2;

	printf("testing line between centre of %d and %d and %d\n",a,b,c);
	if (((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ)) > 100)
	{
		return;
	}


	double jumpScaler = res / ((highX - lowX)*(highX - lowX) + (highY - lowY)*(highY - lowY) + (highZ - lowZ)*(highZ - lowZ));

	double dx = (highX - lowX) * jumpScaler;
	double dy = (highY - lowY) * jumpScaler;
	double dz = (highZ - lowZ) * jumpScaler;
	int reps = 1 / jumpScaler;
	bool sucsess = false;

	printf("%d\n",reps);
	for (size_t i = 0; i < reps; i++)
	{
		int k = i;
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

void runAll(double res, double cutoff,std::string outputfile,int size, wfnData* inputFile,bool makeCube)
{



	//set up threadpool
	boost::asio::io_service ioService;
	boost::thread_group threadpool;
	std::auto_ptr<boost::asio::io_service::work> work(new boost::asio::io_service::work(ioService));

	int numOfThreads = std::thread::hardware_concurrency();
	std::cout << numOfThreads << std::endl;

	for(int i = 0; i< numOfThreads * 2; i++)
		threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));


	for (size_t i = 0; i < (*inputFile).nuc; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			pdrawArgs *lineData;
			lineData = new pdrawArgs(i, j, res, cutoff, outputfile, size, inputFile, makeCube);
			ioService.post(boost::bind(pDrawline, (void *)lineData));

		}
	}

	work.reset();
	threadpool.join_all();
	ioService.stop();

}

std::vector<std::string> readFileLines(const char* filename)
{
	std::vector<std::string> file;
	std::ifstream input(filename);
	std::string line;
	int i = 0;
	while (getline(input, line)){
		file.push_back(line);
	}
	return file;
}
void useInputFile(char* filename)
{

	std::fstream inputFileTest(filename);
	if(!inputFileTest)
	{
		std::cout << "input file not found" << std::endl;
		return;
	}
	std::vector<std::string> lines;
	lines = readFileLines(filename);
	int lineNum = lines.size();
	if(lineNum == 0)
	{
		std::cout << "the input file needs text" <<  std::endl;
		printf("bonder h for help\n");
		return;
	}

	if(lineNum == 1)
	{
		std::cout << "please select option in the input file and the name of the wfn file";
		return;
	}

	wfnData *inputFile = 0;
	try
	{
		inputFile = init(lines[1]);

	}
	catch (const std::invalid_argument& ia) 
	{
		std::cout << "error in parssing wavefunction data, if you have more than 100 atoms 'bonder fixwfn' must be run" << std::endl;
		return;
	}

	if(lines[0] == "p")
	{
		if (lineNum != 9)
		{
			std::cout << "error in parsing input file\npoint file format is:\np\nwfn file\nx\ny\nz\nrdg cutoff\nres\noutput fiel name\noutput cube file"<< std::endl;
			return;
		}
		bool sucsess;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		try
		{
			analize.setUpAnalysisBatch( std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4]), std::stod(lines[5]),batch);
			printf("%f \n", (*batch).RDG(std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4])));
			analize.anilizePoint(0, 0, 0, 0, SIZE, SIZE, std::stod(lines[6]), &sucsess, inputFile, lines[7], batch,lines[8] != "true");
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}

		if (sucsess)
		{
			printf("point given is in region\n");
		}
		else
		{
			printf("point given is not in region\n");
		}
		return;

	}

	if (lines[0] == "l")
	{
		if (lineNum != 8)
		{
			std::cout << "error in parsing input file\nline file format is:\nl\nwfn file\natom1\natom2\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}   
		try
		{
			drawline(std::stoi(lines[2]), std::stoi(lines[3]), std::stod(lines[4]), std::stod(lines[5]), lines[6], SIZE, inputFile, lines[7] != "true");
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}
	//letter file 1 2 res cutoff
	if (lines[0] == "t")
	{
		if (lineNum != 9)
		{
			std::cout << "error in parsing input file\ntriangle file format is:\nt\nwfn file\natom1\natom2\natom3\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}

		try
		{
			drawtrig(std::stoi(lines[2]), std::stoi(lines[3]),std::stoi(lines[4]), std::stod(lines[5]), std::stod(lines[6]), lines[7], SIZE, inputFile, lines[8] != "true");
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}

	if (lines[0] == "q")
	{
		if (lineNum != 10)
		{
			std::cout << "error in parsing input file\ntriangle file format is:\nt\nwfn file\natom1\natom2\natom3\natom4\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}

		try
		{
			drawquad(std::stoi(lines[2]), std::stoi(lines[3]),std::stoi(lines[4]),std::stoi(lines[5]), std::stod(lines[6]), std::stod(lines[7]), lines[8], SIZE, inputFile, lines[9] != "true");
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;

	}
	//letter file res cutoff output
	if (lines[0] == "a")
	{
		if (lineNum != 6)
		{
			std::cout << "error in parsing input file\nall bonds file format is:\na\nwfn file\nrdg cutoff\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}


		try
		{
			runAll(std::stod(lines[2]), std::stod(lines[3]), lines[4], SIZE, inputFile, lines[5] != "true");
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		return;
	}

	//letter file minx miny minz maxx maxy maxz res outputFile
	if (lines[0] == "g")
	{
		if (lineNum != 10)
		{
			std::cout << "error in parsing input file\ngrid file format is:\ng\nwfn file\nlow x\nlow y\n low z\n high x\n high y \n high z\nres\noutput file name\noutput cube file"<< std::endl;
			return;
		}

		analysisBatch* batch = new analysisBatch(*inputFile);
		try
		{

			outputCube(std::stod(lines[2]), std::stod(lines[3]), std::stod(lines[4]), std::stod(lines[5]), std::stod(lines[6]), std::stod(lines[7]), std::stod(lines[8]), lines[9], *inputFile, 1.0, batch, true);
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return;
		}
		printf("done");
		return;
	}

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

	if (argv[1][0] == 'f')
	{
		useInputFile(argv[2]);
		return 0;
	}

	wfnData *inputFile = 0;
	if (argc != 2)
	{
		try
		{
			inputFile = init(argv[2]);

		}
		catch (const std::invalid_argument& ia) 
		{
			std::cout << "error in parssing wavefunction data, if you have more than 100 atoms 'bonder fixwfn' must be run" << std::endl;
			return 1;
		}
	}

	std::cout << "data read" << std::endl;
	//letter file x y z res cutoff
	if (argv[1][0] == 'p')
	{
		if (argc != 9 && argc != 10 )
		{
			printf("arguments are bonder p inputFile x y z res cutoff outputFile\n");
			return 0;
		}
		bool sucsess;
		analysisBatch* batch = new analysisBatch(*inputFile);
		analysis analize = analysis();
		try
		{
			analize.setUpAnalysisBatch( std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]),batch);
			printf("%f \n", (*batch).RDG(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5])));
			if (argc == 10)
			{
				analize.anilizePoint(0, 0, 0, 0, SIZE, SIZE, std::stod(argv[7]), &sucsess, inputFile, argv[8], batch, !strcmp(argv[9], "true"));
			}
			else
			{
				analize.anilizePoint(0, 0, 0, 0, SIZE, SIZE, std::stod(argv[7]), &sucsess, inputFile, argv[8], batch, true);
			}
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
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
			printf("arguments are bonder l inputfile atom1 atom2 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 8)
				drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], SIZE, inputFile,true);
			else
				drawline(std::stoi(argv[3]), std::stoi(argv[4]), std::stod(argv[5]), std::stod(argv[6]), argv[7], SIZE, inputFile, !strcmp(argv[8], "true"));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;

	}

	//letter file 1 2 res cutoff
	if (argv[1][0] == 't')
	{
		if (!(argc == 8 || argc == 9))
		{
			printf("arguments are bonder t inputfile atom1 atom2 atom3 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 9)
				drawtrig(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), argv[8], SIZE, inputFile,true);
			else
				drawtrig(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), argv[8], SIZE, inputFile, !strcmp(argv[9], "true"));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
		return 0;

	}
	
	if (argv[1][0] == 'q')
	{
		if (!(argc == 9 || argc == 10))
		{
			printf("arguments are bonder t inputfile atom1 atom2 atom3 res cutoff outputfile\n");
			return 0;
		}
		try
		{
			if (argc == 9)
				drawquad(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]), std::stod(argv[7]), std::stod(argv[8]), argv[1], SIZE, inputFile,true);
			else
				drawquad(std::stoi(argv[3]), std::stoi(argv[4]),std::stoi(argv[5]),std::stoi(argv[6]), std::stod(argv[7]), std::stod(argv[8]), argv[1], SIZE, inputFile, !strcmp(argv[9], "true"));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
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

		try
		{
			if (argc == 6)
				runAll(std::stod(argv[3]), std::stod(argv[4]), argv[5], SIZE, inputFile,true);
			else
				runAll(std::stod(argv[3]), std::stod(argv[4]), argv[5], SIZE, inputFile, !strcmp(argv[6], "true"));
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
		}
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
		try
		{
			if (argc == 12)
			{
				outputCube(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), argv[10], *inputFile, 1.0, batch, !strcmp(argv[11], "true"));
			}
			else
			{

				outputCube(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), argv[10], *inputFile, 1.0, batch, true);
			}
		}
		catch(const std::invalid_argument& ia)
		{
			std::cout << "error in arguments" << std::endl;
			return 1;
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

