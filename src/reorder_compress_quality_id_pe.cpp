#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include <chrono> //for timing	
#include <omp.h>
#include "sam_block.h"
#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

uint8_t paired_id_code;
std::string preserve_order;
uint32_t numreads, numreads_by_2;
int readlen, num_thr, quantization_level;
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
void encode(FILE *fout, struct qv_options_t *opts, uint32_t readlen, uint32_t numreads, char *quality_array, std::string &infile_order, uint64_t startpos);


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
	num_thr = atoi(argv[3]);
	preserve_order = std::string(argv[4]);	
//	quantization_level = atoi(argv[3]);

	omp_set_num_threads(num_thr);
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.read((char*)&paired_id_code,sizeof(uint8_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;
	generate_order();
	auto start_quality = std::chrono::steady_clock::now();
	reorder_quality();
	auto end_quality = std::chrono::steady_clock::now();
	auto diff_quality = std::chrono::duration_cast<std::chrono::duration<double>>(end_quality-start_quality);
	std::cout << "\nQuality compression total time: " << diff_quality.count() << " s\n";
	auto start_id = std::chrono::steady_clock::now();
	reorder_id();
	auto end_id = std::chrono::steady_clock::now();
	auto diff_id = std::chrono::duration_cast<std::chrono::duration<double>>(end_id-start_id);
	std::cout << "\nID compression total time: " << diff_id.count() << " s\n";
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
	char *quality = new char [(uint64_t)numreads_by_2*(readlen+1)];
	std::string infile_quality[2] = {infile_quality_1,infile_quality_2};		
	for(int k = 0; k < 2; k++)
	{
		std::ifstream f_in(infile_quality[k]);

		for (uint64_t i = 0; i < numreads_by_2; i++)
			f_in.getline((quality+i*(readlen+1)),readlen+1);
		f_in.close();
		#pragma omp parallel	
		{
		int tid = omp_get_thread_num();
		uint64_t start = uint64_t(tid)*(numreads_by_2/omp_get_num_threads());
		uint32_t numreads_thr = numreads_by_2/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			numreads_thr = numreads_by_2-numreads_thr*(omp_get_num_threads()-1);
//		std::ifstream f_order(outfile_order,std::ios::binary);
//		f_order.seekg(start*sizeof(uint32_t));
		struct qv_options_t opts;
		opts.verbose = 1;
		opts.stats = 0;
		opts.ratio = 8.0;
		opts.clusters = 1;
		opts.uncompressed = 0;
		opts.distortion = DISTORTION_MSE;
		opts.cluster_threshold = 4;
		std::string output_name_str = (basedir+"/compressed_quality_"+std::to_string(k+1)+".bin."+std::to_string(tid));
		const char *output_name = output_name_str.c_str();
		FILE *fout;
		fout = fopen(output_name, "wb");
		opts.mode = MODE_FIXED;
		encode(fout, &opts, readlen, numreads_thr, quality, outfile_order,start*sizeof(uint32_t));	
//		f_order.close();
		}
	}
	delete[] quality;
	return;
}

void reorder_id()
{
	std::string *id = new std::string [numreads_by_2];
	std::string infile_id[2] = {infile_id_1,infile_id_2};		
	for(int k = 0; k < 2; k++)
	{
		if(paired_id_code !=0 && k == 1)
			break;	
		std::ifstream f_in(infile_id[k]);

		for (uint64_t i = 0; i < numreads_by_2; i++)
			std::getline(f_in,id[i]);
		f_in.close();
		//facing issues with parallel id compression/decompression
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		uint64_t start = uint64_t(tid)*(numreads_by_2/omp_get_num_threads());
		uint32_t numreads_thr = numreads_by_2/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			numreads_thr = numreads_by_2-numreads_thr*(omp_get_num_threads()-1);	
		std::ifstream f_order(outfile_order,std::ios::binary);
		f_order.seekg(start*sizeof(uint32_t));
		std::string outfile_compressed_id_str = (basedir+"/compressed_id_"+std::to_string(k+1)+".bin."+std::to_string(tid));
		const char *outfile_compressed_id = outfile_compressed_id_str.c_str();
		struct compressor_info_t comp_info;
		comp_info.id_array = id;
		comp_info.f_order = &f_order;
		comp_info.numreads = numreads_thr;
		comp_info.mode = COMPRESSION;
		comp_info.fcomp = fopen(outfile_compressed_id, "w");
		compress((void *)&comp_info);
		fclose(comp_info.fcomp);
		f_order.close();	
		}	

	}
	delete[] id;
	return;
}
