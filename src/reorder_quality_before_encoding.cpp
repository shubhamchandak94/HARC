//Write reordered qualities after the reordering stage. 
//Output written to two files:
//1. matched.quality for qualities of non-singleton reads after reordering
//2. singleton.quality for qualities of singleton reads

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"

uint32_t numreads_matched, numreads_singleton;

std::string outfile;
std::string outfile_s;
std::string infile;
std::string infile_order;
std::string infile_order_s;

void reorder_quality();

void getDataParams();//populate numreads_matched and numreads_singleton

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/matched.quality";
	outfile_s = basedir + "/output/singleton.quality";
	infile = basedir + "/output/input_clean.quality";
	infile_order = basedir + "/output/read_order.bin";
	infile_order_s = basedir + "/output/read_order.bin.singleton";
	getDataParams();
	reorder_quality();
	return 0;
}

void reorder_quality()
{
	std::ofstream f(outfile);
	std::ofstream f_s(outfile_s);
	std::ifstream f_in(infile);
	std::ifstream f_order(infile_order,std::ios::binary);
	std::ifstream f_order_s(infile_order_s,std::ios::binary);

	uint32_t numreads_total = numreads_matched+numreads_singleton;
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads_total];
	for (uint32_t i = 0; i < numreads_matched; i++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order.close();
	
	for (uint32_t i = 0; i < numreads_singleton; i++)
	{
		f_order_s.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i + numreads_matched;
	}
	f_order_s.close();

	//writing qualities of matched reads first
	uint32_t max_bin_size = numreads_matched/4;
	char s[readlen+1];
	s[readlen] = '\0';
	if(max_bin_size != 0)
	{
		for (uint32_t i = 0; i <= numreads_matched/max_bin_size; i++)
		{
			auto numreads_bin = max_bin_size;
			if (i == numreads_matched/max_bin_size)
				numreads_bin = numreads_matched%max_bin_size;
			uint32_t *index_array = new uint32_t [numreads_bin];
			char(*quality_bin)[readlen+1] = new char [numreads_bin][readlen+1];	
			uint32_t pos = 0;
			for(uint32_t j = 0; j < numreads_total; j++)
			{
				order = reverse_index[j];
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
		}
	}	
	//writing qualities of singleton reads
	uint32_t *index_array = new uint32_t [numreads_singleton];
	char(*quality_bin)[readlen+1] = new char [numreads_singleton][readlen+1];	
	uint32_t pos = 0;
	for(uint32_t j = 0; j < numreads_total; j++)
	{
		order = reverse_index[j];
		if (order >= numreads_matched)
		{
			index_array[order-numreads_matched] = pos;
			f_in.seekg(uint64_t(j)*(readlen+1), f_in.beg);
			f_in.getline(quality_bin[pos],readlen+1);
			pos++;
		}
	}
	for(uint32_t j = 0; j < numreads_singleton; j++)
	{
		f_s << quality_bin[index_array[j]] << "\n";
	}
	delete[] index_array;
	delete[] quality_bin;
	f_in.close();
	delete[] reverse_index;
	f.close();
	f_s.close();
	return;
}

void getDataParams()
{
	numreads_matched = 0;
	numreads_singleton = 0;
	std::ifstream f_order(infile_order,std::ios::binary);
	uint32_t order;
	f_order.read((char*)&order,sizeof(uint32_t));
	while(!f_order.eof())
	{
		numreads_matched++;
		f_order.read((char*)&order,sizeof(uint32_t));
	}
	f_order.close();
	
	f_order.open(infile_order_s,std::ios::binary);
	f_order.read((char*)&order,sizeof(uint32_t));
	while(!f_order.eof())
	{
		numreads_singleton++;
		f_order.read((char*)&order,sizeof(uint32_t));
	}
	f_order.close();
	return;
}
