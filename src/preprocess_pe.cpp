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
std::string outfile_meta;
std::string preserve_quality;

int max_readlen=-1, quantization_level;
//quantization level:
//0:lossless
//1:Illumina 8 binning
//2:Binary thresholding

char quantization_table[128];

void generate_quantization_table();

int preprocess();

int quantize(std::string &line, int &readlen);

uint8_t find_id_pattern(std::string &id_1,std::string &id_2);

bool check_id_pattern(std::string &id_1,std::string &id_2, uint8_t paired_id_code);

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[3]);
	infile[0] = std::string(argv[1]);
	infile[1] = std::string(argv[2]);
	preserve_quality = std::string(argv[4]);
//	quantization_level = atoi(argv[6]);
	outfileclean = basedir + "/input_clean.dna";
	outfileN = basedir + "/input_N.dna";
	outfileorderN = basedir + "/read_order_N.bin";
	outfilequality[0] = basedir + "/input_1.quality";
	outfilequality[1] = basedir + "/input_2.quality";
	outfileid[0] = basedir + "/input_1.id";
	outfileid[1] = basedir + "/input_2.id";
	outfilenumreads = basedir + "/numreads.bin";
	outfile_meta = basedir + "/read_meta.txt";
	generate_quantization_table();
	int status = preprocess();
	if(status != 0)
		return -1;
	std::cout << "Preprocessing Done!\n";
	return 0;
}

int preprocess()
{
	std::string line, id_1;
	std::ofstream f_clean(outfileclean);
	std::ofstream f_N(outfileN);
	std::ofstream f_order_N(outfileorderN,std::ios::binary);

	uint64_t total_reads[2] = {0,0};
	uint64_t readnum = 0, num_clean = 0;
	uint8_t paired_id_code = 0;		
	bool paired_id_match = false;
	int current_readlen;
	//code 0: no pattern found
	//code 1: */1 and */2 where * are same in both
	//code 2: * and * where * are same in both
	//code 3: * 1:# and * 2:# where * and # are common to both and * contains no space (used in new versions)
	for(int j = 0; j < 2; j++)
	{
		std::ifstream myfile(infile[j], std::ifstream::in);

		std::ofstream f_quality;
		std::ofstream f_id;
		std::ifstream fin_id_1;

		if(preserve_quality == "True")
		{
			f_quality.open(outfilequality[j]);
			f_id.open(outfileid[j]);
			if(j==1)
			{
				fin_id_1.open(outfileid[0]);
				//check first ids to detect patterns
				std::string id_2;
				std::getline(fin_id_1, id_1);
				std::getline(myfile, id_2);
				paired_id_code = find_id_pattern(id_1,id_2);
				if(paired_id_code != 0)
					paired_id_match = true;
				myfile.close();
				myfile.open(infile[1]);
				fin_id_1.close();
				fin_id_1.open(outfileid[0]);	
			}
		}	
		int i = 0;
		bool flag_N = false;
		while(std::getline(myfile, line))
		{
			switch(i)
			{
				case 0:	if(preserve_quality == "True")
					{
						f_id << line << "\n";
						if(j==1 && paired_id_match)
						{
							std::getline(fin_id_1,id_1);
							if(fin_id_1.eof())
								paired_id_match = false;
							else	
								paired_id_match = check_id_pattern(id_1,line,paired_id_code);
						}
					}
					break;
				case 1: //f << line << "\n";
					current_readlen = line.length();
					if(current_readlen >= 256)
					{
						std::cout << "Read length cannot exceed 255. Read with length "<<current_readlen << " found\n";
						return -1;
					}
					if(current_readlen > max_readlen)
						max_readlen = current_readlen;
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
					{
						if(line.length() != current_readlen)
						{
							std::cout << "Quality length does not match read length: "<< current_readlen << " and " << line.length() << " found.\n";
							return -1;
						}
				/*		if(quantize(line,current_readlen) == -1)
						{
							std::cout << "Invalid quality value found [outside 0 - 41].\n";
						}	
				*/		f_quality << line << "\n";
					}
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
		if(paired_id_match == true)
			f_numreads.write((char*)&paired_id_code, sizeof(uint8_t));
		else
		{
			paired_id_code = 0;
			f_numreads.write((char*)&paired_id_code, sizeof(uint8_t));
		}	
		std::cout << "Max Read length: " << max_readlen << "\n";
		std::cout << "Total number of reads: " << readnum <<"\n";
		std::cout << "Total number of reads without N: " << num_clean <<"\n";
		std::cout << "Paired id match code: " << (int)paired_id_code << "\n";
		f_numreads.close();
		std::ofstream f_meta(outfile_meta);
		f_meta << max_readlen << "\n";
		f_meta.close();
	}	
	return 0;	
}

uint8_t find_id_pattern(std::string &id_1,std::string &id_2)
{
	if(id_1.length() != id_2.length())
		return 0;
	if(id_1 == id_2)
		return 2;
	int len = id_1.length();
	if(id_1[len-1] == '1' && id_2[len-1] == '2')
	{
		//compare rest
		int i;
		for(i = 0; i < len-1; i++)
			if(id_1[i] != id_2[i])
				break;
		if(i == len-1)
			return 1;
	}
	int i = 0;
	for(i = 0; i < len; i++)
	{
		if(id_1[i] != id_2[i])
			break;
		if(id_1[i] == ' ')
		{
			if(i < len-1 && id_1[i+1] == '1' && id_2[i+1] == '2')
				i++;
			else
				break;
		}
	}
	if(i == len)
		return 3;
	return 0;
}

bool check_id_pattern(std::string &id_1,std::string &id_2, uint8_t paired_id_code)
{
	if(id_1.length() != id_2.length())
		return false;
	int len = id_1.length();
	switch(paired_id_code)
	{
		case 1: if(id_1[len-1] == '1' && id_2[len-1] == '2')
			{
				//compare rest
				int i;
				for(i = 0; i < len-1; i++)
					if(id_1[i] != id_2[i])
						break;
				if(i == len-1)
					return true;
			}
			break;
		case 2: if(id_1 == id_2)
				return true;
			break;
		case 3:	int i = 0;
			for(i = 0; i < len; i++)
			{
				if(id_1[i] != id_2[i])
					break;
				if(id_1[i] == ' ')
				{
					if(i < len-1 && id_1[i+1] == '1' && id_2[i+1] == '2')
						i++;
					else
						break;
				}
			}
			if(i == len)
				return true;	  
			break;
	}
	return false;
}

int quantize(std::string &line, int &readlen)
{
	if(quantization_level == 0)
	{
		for(int i = 0; i < readlen; i++)
			if(line[i] > 33+41 || line[i] < 33)//TODO update later
				return -1;
	}
	else
		for(int i = 0; i < readlen; i++)
			line[i] = quantization_table[line[i]];
	return 0;
}

void generate_quantization_table()
{
	switch(quantization_level)
	{
		case 1: for(uint8_t i = 0; i <= 33+9; i++)
				quantization_table[i] = 33+6;
			for(uint8_t i = 33+10; i <= 33+19; i++)
				quantization_table[i] = 33+15;
			for(uint8_t i = 33+20; i <= 33+24; i++)
				quantization_table[i] = 33+22;
			for(uint8_t i = 33+25; i <= 33+29; i++)
				quantization_table[i] = 33+27;
			for(uint8_t i = 33+30; i <= 33+34; i++)
				quantization_table[i] = 33+33;
			for(uint8_t i = 33+35; i <= 33+39; i++)
				quantization_table[i] = 33+37;
			for(uint8_t i = 33+40; i <= 127; i++)
				quantization_table[i] = 33+40;
			break;
		case 2: for(uint8_t i = 0; i <= 33+20; i++)
				quantization_table[i] = 33+6;
			for(uint8_t i = 33+21; i <= 127; i++)
				quantization_table[i] = 33+40;
			break;
	}
	return;		
}
