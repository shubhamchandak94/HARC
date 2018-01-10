#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdio>

std::string infile_order_N;
std::string infilenumreads;
std::string infile_order_paired;
std::string infile_paired_flag_first;
std::string infile_paired_flag_N;

std::string outfile_order;

uint32_t numreads_total, numreads_clean, numreads_N, numreads_by_2; 

void pe_decode();
//create read_order.bin from the paired order and flag files

void unpackbits();
//unpack flag files into 1 byte per flag

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infilenumreads = basedir + "/output/numreads.bin";
	outfile_order = basedir + "/output/read_order.bin";
	infile_order_N = basedir + "/output/read_order_N.bin";
	infile_order_paired = basedir + "/output/read_order_paired.bin";
	infile_paired_flag_first = basedir + "/output/read_paired_flag_first.bin";
	infile_paired_flag_N = basedir + "/output/read_paired_flag_N.bin";	
	
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads_clean,sizeof(uint32_t));
	f_numreads.read((char*)&numreads_total,sizeof(uint32_t));
	f_numreads.close();
	numreads_by_2 = numreads_total/2;
	numreads_N = numreads_total - numreads_clean;	
	
	unpackbits();
	pe_decode();
	return 0;
}

void pe_decode()
{
	std::ifstream in_order_paired(infile_order_paired,std::ios::binary);
	std::ifstream in_paired_flag_first(infile_paired_flag_first);
	std::ifstream in_paired_flag_N(infile_paired_flag_N);
	std::ifstream fin_order_N(infile_order_N,std::ios::binary);
	//first mark locations of pairs containing 'N' reads
	//these locations can't be used later for assigning positions 
	std::vector<bool> pairs_with_N(numreads_by_2,false);
	char flag_N;
	uint32_t order_paired;
	for(uint32_t i = 0; i < numreads_N; i++)
	{
		fin_order_N.read((char*)&order_paired, sizeof(uint32_t));
		pairs_with_N[order_paired%numreads_by_2] = true;
	}
	fin_order_N.close();
	
	//now create index from reordered position to new intended position (numreads_clean --> numreads_total)
	//size numreads_total because we will reuse it later
	uint32_t *read_order = new uint32_t[numreads_total];
	std::fill(read_order, read_order+numreads_total, numreads_total);
	char flag_first;
	uint32_t current_pair = 0;//to track position of available pair
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		while(pairs_with_N[current_pair]==true)
			current_pair++;
		in_paired_flag_first.get(flag_first);
		in_paired_flag_N.get(flag_N);
		if(flag_N == '1')
		{
			in_order_paired.read((char*)&order_paired, sizeof(uint32_t));
			if(flag_first == '1')
				read_order[i] = order_paired;
			else
				read_order[i] = order_paired + numreads_by_2;
		}
		else
		{
			if(read_order[i] != numreads_total)
				continue;//already done by pair
			in_order_paired.read((char*)&order_paired, sizeof(uint32_t));
			if(flag_first == '1')
			{
				read_order[i] = current_pair;
				read_order[i+order_paired] = current_pair + numreads_by_2;
				current_pair++;
			}
			else
			{
				read_order[i] = current_pair + numreads_by_2;
				read_order[i+order_paired] = current_pair;
				current_pair++;			
			}	
		}
	}
	//write read_order array to disk to make space for inverse_order array
	std::ofstream f_temp(outfile_order,std::ios::binary);
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		f_temp.write((char*)&read_order[i], sizeof(uint32_t));
	}
	f_temp.close();
	//now generate inverse index : numreads_total --> numreads_clean
	//ultimately we'll collapse this inverse index to numreads_clean --> numreads_clean
	std::ifstream fin_temp(outfile_order,std::ios::binary);
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		read_order[order_paired] = i;
	}
	fin_temp.close();
		
	//write read_order array to disk to make space for final array
	f_temp.open(outfile_order,std::ios::binary);
	for(uint32_t i = 0; i < numreads_total; i++)
	{
		f_temp.write((char*)&read_order[i], sizeof(uint32_t));
	}
	f_temp.close();
	
	//now use inverse index stored in file along with read_order_N.bin
	//to remove gaps and get a map from numreads_clean --> numreads_clean
	//basically ignore inverse index for N reads

	fin_order_N.open(infile_order_N,std::ios::binary);
	fin_temp.open(outfile_order,std::ios::binary);
	uint32_t pos_in_clean = 0, next_N_read;
	fin_order_N.read((char*)&next_N_read,sizeof(uint32_t));
	for(uint32_t i = 0; i < numreads_total; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		if(i == next_N_read && !fin_order_N.eof())
		{
			fin_order_N.read((char*)&next_N_read,sizeof(uint32_t));
			continue;
		}
		else
		{
			read_order[order_paired] = pos_in_clean;
			pos_in_clean++;
		}			
	}
	fin_order_N.close();
	fin_temp.close();
	
	//finally, write read_order to output file
	std::ofstream fout_order(outfile_order,std::ios::binary);
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		fout_order.write((char*)&read_order[i], sizeof(uint32_t));
	}
	fout_order.close();
	
	return;
}	

void unpackbits()
{
	//flag_first
	std::ifstream in_flag_first(infile_paired_flag_first,std::ios::binary);
	std::ofstream f_flag_first(infile_paired_flag_first+".tmp");
	std::ifstream in_flag_first_tail(infile_paired_flag_first+".tail");
	char inttochar[2];
	inttochar[0] = '0';
	inttochar[1] = '1';
	
	uint8_t packedchar;
	in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	while(!in_flag_first.eof())
	{	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	}
	in_flag_first.close();
	f_flag_first << in_flag_first_tail.rdbuf();
	in_flag_first_tail.close();
	f_flag_first.close();
	infile_paired_flag_first = infile_paired_flag_first + ".tmp";
	
	//flag_N
	std::ifstream in_flag_N(infile_paired_flag_N,std::ios::binary);
	std::ofstream f_flag_N(infile_paired_flag_N+".tmp");
	std::ifstream in_flag_N_tail(infile_paired_flag_N+".tail");
	
	in_flag_N.read((char*)&packedchar,sizeof(uint8_t));
	while(!in_flag_N.eof())
	{	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_N << inttochar[packedchar%2];
		packedchar/=2; 	
		in_flag_N.read((char*)&packedchar,sizeof(uint8_t));
	}
	in_flag_N.close();
	f_flag_N << in_flag_N_tail.rdbuf();
	in_flag_N_tail.close();
	f_flag_N.close();
	infile_paired_flag_N = infile_paired_flag_N + ".tmp";
	return;
}
