#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>

uint32_t numreads;
int readlen;

std::string outfile_quality;
std::string infile_quality;
std::string outfile_id;
std::string infile_id;
std::string infile_order;
std::string infilenumreads;

void reorder_quality();
void reorder_id();

void getDataParams();//populate numreads and numreads_N

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile_quality = basedir + "/output.quality";
	infile_quality = basedir + "/input.quality";
	outfile_id = basedir + "/output.id";
	infile_id = basedir + "/input.id";
	infile_order = basedir + "/read_order.bin";
	infilenumreads = basedir + "/numreads.bin";

	readlen = atoi(argv[2]);
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
	reorder_quality();
	reorder_id();
	return 0;
}

void reorder_quality()
{
	std::ofstream f(outfile_quality);
	std::ifstream f_in(infile_quality);
	std::ifstream f_order(infile_order,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads];
	for (uint32_t i = 0; i < numreads; i++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order.close();
	uint32_t max_bin_size = numreads/4;
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		char *quality_bin = new char [numreads_bin*(readlen+1)];	
		uint32_t pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			order = reverse_index[j];
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f_in.seekg(uint64_t(j)*(readlen+1), f_in.beg);
				f_in.getline((quality_bin+pos*(readlen+1)),readlen+1);
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			f << (quality_bin+index_array[j]*(readlen+1)) << "\n";
		}
		delete[] index_array;
		delete[] quality_bin;
	}
	f_in.close();
	delete[] reverse_index;
	f.close();
	return;
}

void reorder_id()
{
	std::ofstream f(outfile_id);
	std::ifstream f_in(infile_id);
	std::ifstream f_order(infile_order,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads];
	for (uint32_t i = 0; i < numreads; i++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order.close();
	uint32_t max_bin_size = numreads/8;
	std::string s;
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		std::string *id_bin = new std::string [numreads_bin];
		uint32_t pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			order = reverse_index[j];
			std::getline(f_in,s);
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				id_bin[pos] = s;
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			f << id_bin[index_array[j]] << "\n";
		}
		delete[] index_array;
		delete[] id_bin;
	}
	f_in.close();
	delete[] reverse_index;
	f.close();
	return;
}
