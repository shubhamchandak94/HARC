#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include "sam_block.h"
#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

uint8_t paired_id_code;
std::string preserve_order;
uint32_t numreads, numreads_by_2;
int readlen, quantization_level;
//quantization level:
//0:lossless
//1:Illumina 8 binning
//2:Binary thresholding

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
void encode(FILE *fout, struct qv_options_t *opts, uint32_t readlen, uint32_t numreads, char *quality_array, std::ifstream *f_order);


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
	preserve_order = std::string(argv[3]);	
//	quantization_level = atoi(argv[3]);
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.read((char*)&paired_id_code,sizeof(uint8_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;
	generate_order();
	reorder_quality();
	reorder_id();
	return 0;
}

void generate_order()
{
	if(preserve_order == "True") //write fake order information in this case to provide common interface
	{
		std::ofstream fout_order(outfile_order,std::ios::binary);
		for(uint32_t i = 0; i < numreads_by_2; i++)
		{
			fout_order.write((char*)&i,sizeof(uint32_t));
		}
		fout_order.close();
	}
	else
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
	}
	return;
}

void reorder_quality()
{
	char *quality = new char [numreads_by_2*(readlen+1)];
	std::string infile_quality[2] = {infile_quality_1,infile_quality_2};		
//	std::string outfile_quality[2] = {outfile_quality_1,outfile_quality_2};
	for(int k = 0; k < 2; k++)
	{
		std::ifstream f_in(infile_quality[k]);
		std::ifstream f_order(outfile_order,std::ios::binary);
		uint32_t order;
		for (uint64_t i = 0; i < numreads_by_2; i++)
			f_in.getline((quality+i*(readlen+1)),readlen+1);
		f_in.close();
		
		struct qv_options_t opts;
		opts.verbose = 1;
		opts.stats = 0;
		opts.ratio = 8.0;
		opts.clusters = 1;
		opts.uncompressed = 0;
		opts.distortion = DISTORTION_MSE;
		opts.cluster_threshold = 4;
		const char *output_name = (basedir+"/compressed_quality_"+std::to_string(k+1)+".bin.0").c_str();
		FILE *fout = fopen(output_name, "wb");
		opts.mode = MODE_FIXED;
		encode(fout, &opts, readlen, numreads_by_2, quality, &f_order);	
		f_order.close();
	}
	delete[] quality;
	return;
}

void reorder_id()
{
	std::string *id = new std::string [numreads_by_2];
	std::string infile_id[2] = {infile_id_1,infile_id_2};		
//	std::string outfile_id[2] = {outfile_id_1,outfile_id_2};
	for(int k = 0; k < 2; k++)
	{
		if(paired_id_code !=0 && k == 1)
			break;	
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