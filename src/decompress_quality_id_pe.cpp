#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <chrono> //for timing
#include <omp.h>
#include "sam_block.h"
#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

std::string infile_id[2];
std::string outfile_id[2];
std::string infile_quality[2];
std::string outfile_quality[2];
std::string infilenumreads;

int readlen, num_thr, num_thr_e;
uint32_t numreads, numreads_by_2;
uint8_t paired_id_code;

void decompress_id();
void decompress_quality();

void decode(char *input_file, char *output_file, struct qv_options_t *opts);

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile_id[0] = basedir + "/compressed_id_1.bin";
	infile_id[1] = basedir + "/compressed_id_2.bin";
	outfile_id[0] = basedir + "/id_1.txt";
	outfile_id[1] = basedir + "/id_2.txt";
	infile_quality[0] = basedir + "/compressed_quality_1.bin";
	infile_quality[1] = basedir + "/compressed_quality_2.bin";
	outfile_quality[0] = basedir + "/quality_1.txt";
	outfile_quality[1] = basedir + "/quality_2.txt";

	infilenumreads = basedir + "/numreads.bin";
	readlen = atoi(argv[2]);
	num_thr = atoi(argv[3]);
	num_thr_e = atoi(argv[4]);

	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.read((char*)&paired_id_code,sizeof(uint8_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;
	
	omp_set_num_threads(num_thr);
	auto start_quality = std::chrono::steady_clock::now();
	decompress_quality();
	auto end_quality = std::chrono::steady_clock::now();
	auto diff_quality = std::chrono::duration_cast<std::chrono::duration<double>>(end_quality-start_quality);
	std::cout << "\nQuality decompression total time: " << diff_quality.count() << " s\n";
	auto start_id = std::chrono::steady_clock::now();
	decompress_id();
	auto end_id = std::chrono::steady_clock::now();
	auto diff_id = std::chrono::duration_cast<std::chrono::duration<double>>(end_id-start_id);
	std::cout << "\nID decompression total time: " << diff_id.count() << " s\n";
}

void decompress_id()
{
	for(int k = 0; k < 2; k++)
	{
		if(paired_id_code != 0 && k==1)
			break;
/*
		struct compressor_info_t comp_info;
		comp_info.numreads = numreads_by_2;
		comp_info.mode = DECOMPRESSION;
		comp_info.fcomp = fopen((infile_id[k]+".0").c_str(),"r");
		comp_info.f_id = fopen((outfile_id[k]+".0").c_str(),"w");
		decompress((void *)&comp_info);
		fclose(comp_info.fcomp);
		fclose(comp_info.f_id);
*/

//facing issues with parallel id

		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{
			uint32_t numreads_thr = numreads_by_2/num_thr_e;
			if(tid_e == num_thr_e - 1)
				numreads_thr = numreads_by_2-numreads_thr*(num_thr_e-1); 
			struct compressor_info_t comp_info;
			comp_info.numreads = numreads_thr;
			comp_info.mode = DECOMPRESSION;
			comp_info.fcomp = fopen((infile_id[k]+"."+std::to_string(tid_e)).c_str(),"r");
			comp_info.f_id = fopen((outfile_id[k]+"."+std::to_string(tid_e)).c_str(),"w");
			decompress((void *)&comp_info);
			fclose(comp_info.fcomp);
			fclose(comp_info.f_id);
		}
		}	

	}
	return;
}

void decompress_quality()
{
	for(int k = 0; k < 2; k++)
	{	
		#pragma omp parallel
		{
		int tid = omp_get_thread_num();
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{		
			struct qv_options_t opts;
			opts.verbose = 1;
			std::string input_file_string = (infile_quality[k]+"."+std::to_string(tid_e));
			char *input_file = new char [input_file_string.length()+1];
			strcpy(input_file,input_file_string.c_str());
			std::string output_file_string = (outfile_quality[k]+"."+std::to_string(tid_e));
			char *output_file = new char [output_file_string.length()+1];
			strcpy(output_file,output_file_string.c_str());
			decode(input_file,output_file,&opts);
		}
		}
	}
	return;
}
