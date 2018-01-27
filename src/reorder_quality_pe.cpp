#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>

uint32_t numreads, numreads_by_2;
int readlen;

std::string outfile_quality_1;
std::string outfile_quality_2;
std::string infile_quality_1;
std::string infile_quality_2;
std::string outfile_id_1;
std::string outfile_id_2;
std::string infile_id_1;
std::string infile_id_2;
std::string infile_order;
std::string outfile_order;
std::string infilenumreads;

void generate_order();
//generate reordering information for the two separate files (pairs) from read_order.bin

void reorder_quality(std::string infile_quality, std::string outfile_quality, std::string infile_order, uint32_t numreads);
void reorder_id(std::string infile_id, std::string outfile_id, std::string infile_order, uint32_t numreads);

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile_quality_1 = basedir + "/output_1.quality";
	outfile_quality_2 = basedir + "/output_2.quality";
	infile_quality_1 = basedir + "/input_1.quality";
	infile_quality_2 = basedir + "/input_2.quality";
	outfile_id_1 = basedir + "/output_1.id";
	outfile_id_2 = basedir + "/output_2.id";
	infile_id_1 = basedir + "/input_1.id";
	infile_id_2 = basedir + "/input_2.id";
	infile_order = basedir + "/read_order.bin";
	outfile_order = basedir + "/read_order.bin.tmp";
	infilenumreads = basedir + "/numreads.bin";

	readlen = atoi(argv[2]);
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;
	generate_order();
	reorder_quality(infile_quality_1,outfile_quality_1,outfile_order,numreads_by_2);
	reorder_quality(infile_quality_2,outfile_quality_2,outfile_order,numreads_by_2);
	reorder_id(infile_id_1,outfile_id_1,outfile_order,numreads_by_2);
	reorder_id(infile_id_2,outfile_id_2,outfile_order,numreads_by_2);
	return 0;
}

void generate_order()
{
	std::ifstream fin_order(infile_order,std::ios::binary);
	std::ofstream fout_order(outfile_order,std::ios::binary);
	uint32_t order;
	for(uint32_t i = 0; i < numreads; i++)
	{
		fin_order.read((char*)&order,sizeof(uint32_t));
		if(order < numreads_by_2)
		{
			fout_order.write((char*)&order,sizeof(uint32_t));
		}
	}
	fin_order.close();
	fout_order.close();
	return;
}

void reorder_quality(std::string infile_quality, std::string outfile_quality, std::string infile_order, uint32_t numreads)
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
	uint64_t max_bin_size = numreads/2;
	for (uint64_t i = 0; i <= numreads/max_bin_size; i++)
	{
		uint64_t numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		char *quality_bin = new char [numreads_bin*(readlen+1)];	
		uint64_t pos = 0;
		for(uint64_t j = 0; j < numreads; j++)
		{
			order = reverse_index[j];
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f_in.seekg(j*(readlen+1), f_in.beg);
				f_in.getline((quality_bin+pos*(readlen+1)),readlen+1);
				pos++;
			}
		}
		for(uint64_t j = 0; j < numreads_bin; j++)
		{
			f << (quality_bin+uint64_t(index_array[j])*(readlen+1)) << "\n";
		}

		delete[] index_array;
		delete[] quality_bin;
	}
	f_in.close();
	delete[] reverse_index;
	f.close();
	return;
}

void reorder_id(std::string infile_id, std::string outfile_id, std::string infile_order, uint32_t numreads)
{
	std::ofstream f(outfile_id);
	std::ifstream f_in;
	std::ifstream f_order(infile_order,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads];
	for (uint32_t i = 0; i < numreads; i++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order.close();
	uint32_t max_bin_size = numreads/2;
	std::string s;
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		f_in.open(infile_id);
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
		f_in.close();
		delete[] index_array;
		delete[] id_bin;
	}
	delete[] reverse_index;
	f.close();
	return;
}
