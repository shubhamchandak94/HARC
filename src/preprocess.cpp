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
std::string outfileid;
std::string outfileidN;
std::string outfileidfinal;
std::string outfilenumreads;

std::string preserve_order;
std::string preserve_quality;

int readlen;


int preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[2]);
	infile = std::string(argv[1]);
	preserve_order = std::string(argv[3]);
	preserve_quality = std::string(argv[4]);
	readlen = atoi(argv[5]);
	outfileclean = basedir + "/output/input_clean.dna";
	outfileN = basedir + "/output/input_N.dna";
	outfileorderN = basedir + "/output/read_order_N.bin";
	outfilequality = basedir + "/output/input_clean.quality";
	outfilequalityN = basedir + "/output/input_N.quality";
	outfilequalityfinal = basedir + "/output/output.quality";
	outfileid = basedir + "/output/input_clean.id";
	outfileidN = basedir + "/output/input_N.id";
	outfileidfinal = basedir + "/output/output.id";
	outfilenumreads = basedir + "/output/numreads.bin";
	int status = preprocess();
	if(status != 0)
		return -1;
	std::cout << "Preprocessing Done!\n";
	return 0;
}

int preprocess()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	std::ofstream f_order_N(outfileorderN,std::ios::binary);
	std::ofstream f_quality;
	std::ofstream f_id;
	
	if(preserve_quality == "True")
		if(preserve_order == "False")	
		{
			f_quality.open(outfilequality);
			f_id.open(outfileid);
		}
		else
		{
			f_quality.open(outfilequalityfinal);
			f_id.open(outfileidfinal);
		}
	std::ofstream f_quality_N;
	std::ofstream f_id_N;
	if(preserve_order == "False" && preserve_quality == "True")
	{
		f_quality_N.open(outfilequalityN);
		f_id_N.open(outfileidN);
	}
	int i = 0;
	uint64_t readnum = 0;
	uint64_t num_clean = 0;
	bool flag_N = false;
	while(std::getline(myfile, line))
	{
		switch(i)
		{
			case 0:	if(preserve_quality == "True")
					if(!flag_N || preserve_order == "True")
						f_id << line << "\n";
					else
						f_id_N << line << "\n";
				break;
			case 1: //f << line << "\n";
				if(line.length() != readlen)
				{	
					std::cout << "Read length not fixed. Found two different read lengths: "<< 
						  readlen << " and " << line.length() << "\n";
					return -1;
				}
				if(line.find('N')!=std::string::npos)
				{
					flag_N = true;
					f_N << line << "\n";
					f_order_N.write((char*)&readnum,sizeof(uint32_t));
				}
				else
				{
					num_clean++;
					flag_N = false;
					f_clean << line << "\n";
				}
				break;
			case 2: break;
			case 3: if(preserve_quality == "True")
					if(!flag_N || preserve_order == "True")
						f_quality << line << "\n";
					else
						f_quality_N << line << "\n";
				readnum++;
				break;
		}
		i = (i+1)%4;
	}
	if(readnum > 4294967290)
	{
		std::cout << "Too many reads. HARC supports at most 4294967290 reads\n";
		return -1;
	}
	else
	{
		std::ofstream f_numreads(outfilenumreads,std::ios::binary);
		uint32_t num_clean_32 = num_clean;
		f_numreads.write((char*)&num_clean_32, sizeof(uint32_t));
		f_numreads.write((char*)&((uint32_t)readnum), sizeof(uint32_t));
		std::cout << "Read length: " << readlen << "\n";
		std::cout << "Total number of reads: " << readnum <<"\n";
		std::cout << "Total number of reads without N: " << num_clean <<"\n";
		f_numreads.close();
	}	
	return 0;	
}
