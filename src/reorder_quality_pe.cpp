#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include "sam_block.h"

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
std::string basedir;
void generate_order();
//generate reordering information for the two separate files (pairs) from read_order.bin

void reorder_quality();
void reorder_id();

int main(int argc, char** argv)
{
	basedir = std::string(argv[1]);
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
	reorder_quality();
	reorder_id();
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

void reorder_quality()
{
	char *quality = new char [numreads_by_2*(readlen+1)];
	std::string infile_quality[2] = {infile_quality_1,infile_quality_2};		
	std::string outfile_quality[2] = {outfile_quality_1,outfile_quality_2};
	for(int k = 0; k < 2; k++)
	{
		std::ofstream f(outfile_quality[k]);
		std::ifstream f_in(infile_quality[k]);
		std::ifstream f_order(outfile_order,std::ios::binary);
		uint32_t order;
		for (uint64_t i = 0; i < numreads_by_2; i++)
			f_in.getline((quality+i*(readlen+1)),readlen+1);
		for (uint64_t i = 0; i < numreads_by_2; i++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			f << (quality+uint64_t(order)*(readlen+1)) << "\n";	
		}
		f_in.close();
		f.close();
		f_order.close();
	}
	delete[] quality;
	return;
}

void reorder_id()
{
	std::string *id = new std::string [numreads_by_2];
	std::string infile_id[2] = {infile_id_1,infile_id_2};		
	std::string outfile_id[2] = {outfile_id_1,outfile_id_2};
	for(int k = 0; k < 2; k++)
	{
		std::ifstream f_in(infile_id[k]);
		std::ifstream f_order(outfile_order,std::ios::binary);
		uint32_t order;
		for (uint64_t i = 0; i < numreads_by_2; i++)
			std::getline(f_in,id[i]);
		f_in.close();
	
		const char *outfile_compressed_id = (basedir+"/compressed_id_"+std::to_string(k+1)+".bin.0").c_str();
		struct compressor_info_t comp_info;
		comp_info.id_array = id;
		comp_info.f_order = &f_order;
		comp_info.numreads = numreads_by_2;
		comp_info.mode = COMPRESSION;
		comp_info.fcomp = fopen(outfile_compressed_id, "w");
		compress((void *)&comp_info);
		fclose(comp_info.fcomp);
		f_order.close();		
	}
	delete[] id;
	return;
}
