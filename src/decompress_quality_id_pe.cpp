#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
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
	f_numreads.close();
	numreads_by_2 = numreads/2;
	
	decompress_quality();
	decompress_id();
}

void decompress_id()
{
	for(int k = 0; k < 2; k++)
	{
		struct compressor_info_t comp_info;
		comp_info.numreads = numreads_by_2;
		comp_info.mode = DECOMPRESSION;
		comp_info.fcomp = fopen((infile_id[k]+".0").c_str(), "r");
		comp_info.f_id = fopen((outfile_id[k]+".0").c_str(),"w");
		decompress((void *)&comp_info);
		fclose(comp_info.fcomp);
		fclose(comp_info.f_id);
	}
	return;
}

void decompress_quality()
{
	for(int k = 0; k < 2; k++)
	{
		struct qv_options_t opts;
		opts.verbose = 1;
		std::string input_file_string = (infile_quality[k]+".0");
		char *input_file = new char [input_file_string.length()+1];
		strcpy(input_file,input_file_string.c_str());
		std::string output_file_string = (outfile_quality[k]+".0");
		char *output_file = new char [output_file_string.length()+1];
		strcpy(output_file,output_file_string.c_str());
		decode(input_file,output_file,&opts);
	}
	return;
}
