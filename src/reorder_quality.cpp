#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"

uint32_t numreads, numreads_N;

std::string outfile;
std::string infile;
std::string infile_N;
std::string outfile_id;
std::string infile_id;
std::string infile_id_N;
std::string infile_order;
std::string infile_order_N_pe;//post-encoding

void reorder_quality();
void reorder_quality_N();
void reorder_id();
void reorder_id_N();

void getDataParams();//populate numreads and numreads_N

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.quality";
	infile = basedir + "/output/input_clean.quality";
	infile_N = basedir + "/output/input_N.quality";
	outfile_id = basedir + "/output/output.id";
	infile_id = basedir + "/output/input_clean.id";
	infile_id_N = basedir + "/output/input_N.id";
	infile_order = basedir + "/output/read_order.bin";
	infile_order_N_pe = basedir + "/output/read_order_N_pe.bin";
	getDataParams();
	reorder_quality_N();
	reorder_quality();
	reorder_id_N();
	reorder_id();
	return 0;
}

void reorder_quality()
{
	std::ofstream f(outfile);
	std::ifstream f_in(infile);
	std::ifstream f_N(infile_N);
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
	char s[readlen+1];
	s[readlen] = '\0';
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		char(*quality_bin)[readlen+1] = new char [numreads_bin][readlen+1];	
		uint32_t pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
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
	f_in.close();
	f << f_N.rdbuf();
	delete[] reverse_index;
	f_N.close();
	f.close();
	return;
}

void reorder_quality_N()
{
	std::ifstream f_N(infile_N);
	std::ifstream f_order_N_pe(infile_order_N_pe,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads_N];
	for (uint32_t i = 0; i < numreads_N; i++)
	{
		f_order_N_pe.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order_N_pe.close();
	char s[readlen+1];
	s[readlen] = '\0';
	uint32_t *index_array = new uint32_t [numreads_N];
	char(*quality)[readlen+1] = new char [numreads_N][readlen+1];	
	uint32_t pos = 0;
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		index_array[reverse_index[j]] = j;
		f_N.getline(quality[j],readlen+1);
	}
	f_N.close();
	std::ofstream f(infile_N);
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		f << quality[index_array[j]] << "\n";
	}
	delete[] index_array;
	delete[] quality;
	delete[] reverse_index;

	f.close();
	return;
}

void reorder_id()
{
	std::ofstream f(outfile_id);
	std::ifstream f_in(infile_id);
	std::ifstream f_N(infile_id_N);
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
	f << f_N.rdbuf();
	delete[] reverse_index;
	f_N.close();
	f.close();
	return;
}

void reorder_id_N()
{
	std::ifstream f_N(infile_id_N);
	std::ifstream f_order_N_pe(infile_order_N_pe,std::ios::binary);
	uint32_t order;
	uint32_t *reverse_index = new uint32_t [numreads_N];
	for (uint32_t i = 0; i < numreads_N; i++)
	{
		f_order_N_pe.read((char*)&order,sizeof(uint32_t));
		reverse_index[order] = i;
	}
	f_order_N_pe.close();
	char s[readlen+1];
	s[readlen] = '\0';
	uint32_t *index_array = new uint32_t [numreads_N];
	std::string *id = new std::string [numreads_N];	
	uint32_t pos = 0;
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		index_array[reverse_index[j]] = j;
		std::getline(f_N,id[j]);
	}
	f_N.close();
	std::ofstream f(infile_N);
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		f << id[index_array[j]] << "\n";
	}
	delete[] index_array;
	delete[] id;
	delete[] reverse_index;

	f.close();
	return;
}

void getDataParams()
{
	uint32_t number_of_lines = 0;
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);
	while (std::getline(myfile, line))
		++number_of_lines;
	numreads = number_of_lines;
	myfile.close();
	number_of_lines = 0;
	std::ifstream myfile_N(infile_N, std::ifstream::in);
	while (std::getline(myfile_N, line))
		++number_of_lines;
	numreads_N = number_of_lines;
	myfile_N.close();
}
