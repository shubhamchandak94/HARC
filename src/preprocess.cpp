#include <iostream>
#include <fstream>
#include <string>

std::string infile;
std::string outfileclean;
std::string outfileN;
std::string outfileorderN;
std::string outfilequality;
std::string outfilequalityN;
std::string outfilequalityfinal;

std::string preserve_order;
std::string preserve_quality;

void preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[2]);
	
	infile = std::string(argv[1]);
	preserve_order = std::string(argv[3]);
	preserve_quality = std::string(argv[4]);
	outfileclean = basedir + "/output/input_clean.dna";
	outfileN = basedir + "/output/input_N.dna";
	outfileorderN = basedir + "/output/read_order_N.bin";
	outfilequality = basedir + "/output/input_clean.quality";
	outfilequalityN = basedir + "/output/input_N.quality";
	outfilequalityfinal = basedir + "/output/output.quality";
	preprocess();
	std::cout << "Preprocessing Done!\n";
	return 0;
}

void preprocess()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	std::ofstream f_order_N(outfileorderN,std::ios::binary);
	std::ofstream f_quality;
	if(preserve_quality == "True")
		if(preserve_order == "False")	
			f_quality.open(outfilequality);
		else
			f_quality.open(outfilequalityfinal);
	std::ofstream f_quality_N;
	if(preserve_order == "False" && preserve_quality == "True")
		f_quality_N.open(outfilequalityN);
	int i = 0;
	uint32_t readnum = 0;
	bool flag_N = false;
	while(std::getline(myfile, line))
	{
		switch(i)
		{
			case 0:	//f_id << line << "\n";
				break;
			case 1: //f << line << "\n";
				if(line.find('N')!=std::string::npos)
				{
					flag_N = true;
					f_N << line << "\n";
					f_order_N.write((char*)&readnum,sizeof(uint32_t));
				}
				else
				{
					flag_N = false;
					f_clean << line << "\n";
				}
				break;
			case 2: break;
			case 3:	//f_quality << line << "\n";
				if(preserve_quality == "True")
					if(!flag_N || preserve_order == "True")
						f_quality << line << "\n";
					else
						f_quality_N << line << "\n";
				readnum++;
				break;
		}
		i = (i+1)%4;
	}
	return;	
}
