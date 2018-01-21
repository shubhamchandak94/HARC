#include <iostream>
#include <fstream>
#include <string>

std::string infile[2];
std::string outfileclean;
std::string outfileN;
std::string outfileorderN;
std::string outfilequality[2];
std::string outfileid[2];
std::string outfilenumreads;

std::string preserve_quality;

int readlen;

int preprocess();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[3]);
	infile[0] = std::string(argv[1]);
	infile[1] = std::string(argv[2]);
	preserve_quality = std::string(argv[4]);
	readlen = atoi(argv[5]);
	outfileclean = basedir + "/input_clean.dna";
	outfileN = basedir + "/input_N.dna";
	outfileorderN = basedir + "/read_order_N.bin";
	outfilequality[0] = basedir + "/input_1.quality";
	outfilequality[1] = basedir + "/input_2.quality";
	outfileid[0] = basedir + "/input_1.id";
	outfileid[1] = basedir + "/input_2.id";
	outfilenumreads = basedir + "/numreads.bin";
	int status = preprocess();
	if(status != 0)
		return -1;
	std::cout << "Preprocessing Done!\n";
	return 0;
}

int preprocess()
{
	std::string line;
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	std::ofstream f_order_N(outfileorderN,std::ios::binary);

	uint64_t total_reads[2] = {0,0};
	uint64_t readnum = 0, num_clean = 0;
		
	for(int j = 0; j < 2; j++)
	{
		std::ifstream myfile(infile[j], std::ifstream::in);

		std::ofstream f_quality;
		std::ofstream f_id;
		
		if(preserve_quality == "True")
		{
			f_quality.open(outfilequality[j]);
			f_id.open(outfileid[j]);
		}	
		int i = 0;
		bool flag_N = false;
		while(std::getline(myfile, line))
		{
			switch(i)
			{
				case 0:	if(preserve_quality == "True")
						f_id << line << "\n";
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
						f_quality << line << "\n";
					readnum++;
					break;
			}
			i = (i+1)%4;
		}
		total_reads[j] = readnum;
	}
	total_reads[1] = total_reads[1] - total_reads[0];	
	if(readnum > 4294967290)
	{
		std::cout << "Too many reads. HARC supports at most 4294967290 reads\n";
		return -1;
	}
	else if(total_reads[1] != total_reads[0])
	{
		std::cout << "Number of reads in the two paired files are not equal\n";
		return -1;
	}
	else
	{
		std::ofstream f_numreads(outfilenumreads,std::ios::binary);
		uint32_t num_clean_32 = num_clean;
		uint32_t readnum_32 = readnum;
		f_numreads.write((char*)&num_clean_32, sizeof(uint32_t));
		f_numreads.write((char*)&readnum_32, sizeof(uint32_t));
		std::cout << "Read length: " << readlen << "\n";
		std::cout << "Total number of reads: " << readnum <<"\n";
		std::cout << "Total number of reads without N: " << num_clean <<"\n";
		f_numreads.close();
	}	
	return 0;	
}
