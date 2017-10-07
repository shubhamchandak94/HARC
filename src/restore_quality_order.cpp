//Restore the order of quality values (order-preserving mode) after encoding.
//Input taken from read_non_N.quality and read_N.quality produced by 
//encoding stage. Output written to output.quality

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include "config.h"

std::string infile_quality;
std::string infile_quality_N;
std::string infile_order;
std::string infile_order_N_pe;
std::string infile_order_N;
std::string infilenumreads;

std::string outfile_quality;

uint32_t numreads = 0;

void separate_N();
void restore_order_N();
void restore_order_non_N();
void merge_N();


int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile_quality = basedir + "/output/read_non_N.quality";
	infile_quality_N = basedir + "/output/read_N.quality";
	infile_order = basedir + "/output/read_order.bin";
	infile_order_N_pe =  basedir + "/output/read_order_N_pe.bin";
	infile_order_N = basedir + "/output/read_order_N.bin";
	
	outfile_quality = basedir + "/output/output.quality";
	infilenumreads = basedir + "/output/numreads.bin";
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));	
	restore_order_N();
	restore_order_non_N();
	merge_N();
	return 0;
}

void restore_order_N()
{	
	std::string line;
	std::ifstream f_N(infile_quality_N);
	uint32_t numreads_N = 0;
	while(std::getline(f_N,line))
		numreads_N++;
	f_N.close();
	f_N.open(infile_quality_N);
		
	std::ifstream f_order(infile_order_N_pe,std::ios::binary);
	uint32_t *index_array = new uint32_t [numreads_N];
	char(*quality_N)[readlen+1] = new char [numreads_N][readlen+1];
	uint32_t order;
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		index_array[order] = j;
		f_N.getline(quality_N[j],readlen+1);
	}
	f_N.close();
	std::ofstream out_N(infile_quality_N+".tmp");
	
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		out_N << quality_N[index_array[j]] << "\n";
	}
	delete[] index_array;
	delete[] quality_N;
	f_order.close();
	out_N.close();
	return;
}

void restore_order_non_N()
{
	uint64_t max_bin_size = numreads/4;
	char s[readlen+1];
	s[readlen] = '\0';
	std::ofstream f(infile_quality+".tmp");
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		std::ifstream f_order(infile_order,std::ios::binary);
		std::ifstream f_in(infile_quality);
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		char(*quality_bin)[readlen+1] = new char [numreads][readlen+1];
		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f_in.seekg(uint64_t(j)*(readlen+1), f_in.beg);
				f_in.getline(quality_bin[pos],readlen+1);
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			f << quality_bin[index_array[j]] << "\n";
		}
		delete[] index_array;
		delete[] quality_bin;
		f_order.close();
		f_in.close();
	}
	f.close();
	return;
}

void merge_N()
{
	std::ofstream f(outfile_quality);
	std::ifstream f_in(infile_quality+".tmp");
	std::ifstream f_in_N(infile_quality_N+".tmp");
	std::ifstream f_order_N(infile_order_N,std::ios::binary);

	std::string line_N,line;
	uint32_t next,i=0;
	f_order_N.read((char*)&next,sizeof(uint32_t));
	while (std::getline(f_in, line))
	{
		while(i==next && !f_order_N.eof())
		{
			std::getline(f_in_N, line_N);
			f << line_N << "\n";
			i++;
			f_order_N.read((char*)&next,sizeof(uint32_t));
		}
		f << line << "\n";
		i++;
	}	
	while(i==next && !f_order_N.eof())
	{
		std::getline(f_in_N, line_N);
		f << line_N << "\n";
		i++;
		f_order_N.read((char*)&next,sizeof(uint32_t));
	}
	f.close();
	f_in.close();
	f_in_N.close();
	f_order_N.close();
	remove((infile_quality+".tmp").c_str());
	remove((infile_quality_N+".tmp").c_str());
	return;
}
