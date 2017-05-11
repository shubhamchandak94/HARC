#include <iostream>
#include <fstream>
#include <string>

std::string infile;
std::string outfile;
std::string outfilequality;
std::string outfileclean;
std::string outfileN;

void preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/input.fastq";
	outfile = basedir + "/input.dna";
	outfilequality = basedir + "/input.quality";
	outfileclean = basedir + "/input_clean.dna";
	outfileN = basedir + "/input_N.dna";
	preprocess();
	std::cout << "Preprocessing Done!\n";
	return 0;
}

void preprocess()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	std::ofstream f(outfile);
	std::ofstream f_quality(outfilequality);
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	
	int i = 0;
	while(std::getline(myfile, line))
	{
		switch(i)
		{
			case 0:	break;
			case 1: f << line << "\n";
				if(line.find('N')!=std::string::npos)
					f_N << line << "\n";
				else
					f_clean << line << "\n";
				break;
			case 2: break;
			case 3:	f_quality << line << "\n";
				break;
		}
		i = (i+1)%4;
	}
	return;	
}
