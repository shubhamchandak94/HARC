#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

std::string infile_order;
std::string infile_order_N;
std::string infilenumreads;

std::string outfile_order_paired;
std::string outfile_paired_flag_first;
std::string outfile_paired_flag_N;

uint32_t numreads_total, numreads_clean, numread_N, numreads_by_2; 

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order, bool* read_flag_N);
//populate arrays:
//read_order = pos in reordered file to pos in original file
//read_inverse_order = pos in original file to pos in reordered file
//read_flag_N = bool array indicating reads which contain N

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order, bool* read_flag_N);
//write to all output files
//order_paired - store relative position of paired read once per pair (store for the read occuring first in the reordered file)
//For a clean read whose pair is an N read, store the actual position and mark in paired_flag_N
//paired_flag_first stores 1 for 1st read of pair and 0 for second read.

void packbits();
//pack flag files into 1 bit per flag

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infilenumreads = basedir + "/output/numreads.bin";
	infile_order = basedir + "/output/read_order.bin";
	infile_order_N = basedir + "/output/read_order_N.bin";
	outfile_order_paired = basedir + "/output/read_order_paired.bin";
	outfile_paired_flag_first = basedir + "/output/read_paired_flag_first.bin";
	outfile_paired_flag_N = basedir + "/output/read_paired_flag_N.bin";	
	
	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.read((char*)&numreads_clean,sizeof(uint32_t));
	f_numreads.read((char*)&numreads_total,sizeof(uint32_t));
	f_numreads.close();
	numread_by_2 = numreads_total/2;
	numreads_N = numreads_total - numreads_N;	
	
	uint32_t *read_order = new uint32_t [numreads_clean];
	uint32_t *read_inverse_order = new uint32_t[numreads_total];
	bool *read_flag_N = new bool[numreads_total]();
	populate_arrays(read_order, read_inverse_order, read_flag_N);
    write_order_paired(read_order, read_inverse_order, read_flag_N);
	packbits();
	
	delete[] read_order;
	delete[] read_inverse_order;
	delete[] read_flag_N;
	return 0;
}

void populate_arrays(uint32_t* read_order, uint32_t* read_inverse_order, bool* read_flag_N)
{
	//read files read_order and read_order_N
	std::ifstream f_order(infile_order, std::ios::binary);
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		f_order.read((char*)&read_order[i], sizeof(uint32_t));
	}
	f_order.close();
	std::ifstream f_order_N(infile_order_N, std::ios::binary);
	uint32_t pos_N;
	for(uint32_t i = 0; i < numreads_N; i++)
	{
		f_order_N.read((char*)&pos_N, sizeof(uint32_t));
		read_flag_N[pos_N] = true;
	}
	f_order_N.close();
	
	//temporarily store cumulative number of N reads in read_inverse_order and 
	//update read_order to contain position in original file rather than clean file.
	uint32_t pos_in_clean = 0, num_N_reads_till_now = 0;
	for(uint32_t i = 0; i < numreads_total; i++)
	{
		if(read_flag_N[i] == true)
			num_N_reads_till_now++;	
		else
			read_inverse_order[j++] = num_N_reads_till_now++;
	}
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		read_order[i] += read_inverse_order[i];
	}

	//now fill read_inverse_order
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		read_inverse_order[read_order[i]] = i;
	}
	return;
}	

void write_order_paired(uint32_t* read_order, uint32_t* read_inverse_order, bool* read_flag_N)	
{
	std::ofstream f_flag_first(outfile_paired_flag_first);
	std::ofstream f_flag_N(outfile_paired_flag_N);
	std::ofstream f_order_paired(outfile_order_paired,std::ios::binary);
	for(uint32_t i = 0; i < numreads_clean; i++)
	{
		if(read_order[i] < numreads_by_2)
		{
			f_flag_first << '1';
			if(read_flag_N[read_order[i] + numreads_by_2])//pair is N read
			{
				f_flag_N << '1';
				f_order_paired.write((char*)&read_order[i],sizeof(uint32_t));
			}
			else if(read_inverse_order[read_order[i] + numreads_by_2] < i)//pair already seen
			{
				f_flag_N << '0';	
			}
			else
			{
				f_flag_N << '0';	
				f_order_paired.write((char*)&(read_inverse_order[read_order[i] + numreads_by_2]-i),sizeof(uint32_t));
			}
		}
		else
		{
			f_flag_first << '0';
			if(read_flag_N[read_order[i] - numreads_by_2])//pair is N read
			{
				f_flag_N << '1';
				f_order_paired.write((char*)&read_order[i],sizeof(uint32_t));
			}
			else if(read_inverse_order[read_order[i] - numreads_by_2] < i)//pair already seen
			{
				f_flag_N << '0';	
			}
			else
			{
				f_flag_N << '0';	
				f_order_paired.write((char*)&(read_inverse_order[read_order[i] - numreads_by_2]-i),sizeof(uint32_t));
			}
		}			
	}
	f_flag_first.close();
	f_flag_N.close();
	f_order_paired.close();
}

void packbits()
{
	//flag_first
	std::ifstream in_flag_first(outfile_paired_flag_first);
	std::ofstream f_flag_first(outfile_paired_flag_first+".tmp",std::ios::binary);
	std::ofstream f_flag_first_tail(outfile_paired_flag_first+".tail");
	
	uint8_t chartoint[128];
	chartoint['0'] = 0;
	chartoint['1'] = 1;
	in_flag_first.close();
	in_flag_first.open(outfile_paired_flag_first);
	char chararray[8];
	uint8_t packedchar;
	for(uint64_t i = 0; i < numreads_clean/8; i++)
	{
		in_flag_first.read(chararray,8);
		
		packedchar = 128*chartoint[chararray[7]]+ 64*chartoint[chararray[6]]
				+ 32*chartoint[chararray[5]]+ 16*chartoint[chararray[4]]
				+  8*chartoint[chararray[3]]+  4*chartoint[chararray[2]]
				+  2*chartoint[chararray[1]]+  1*chartoint[chararray[0]];
		f_flag_first.write((char*)&packedchar,sizeof(uint8_t));
	}
	f_flag_first.close();
	in_flag_first.read(chararray,numreads_clean%8);
	for(int i=0; i<numreads_clean%8;i++)
		f_flag_first_tail << chararray[i];
	f_flag_first_tail.close();
	in_flag_first.close();
	remove((outfile_paired_flag_first).c_str());
	rename((outfile_paired_flag_first+".tmp").c_str(),(outfile_paired_flag_first).c_str());		

	//flag_N
	std::ifstream in_flag_N(outfile_paired_flag_N);
	std::ofstream f_flag_N(outfile_paired_flag_N+".tmp",std::ios::binary);
	std::ofstream f_flag_N_tail(outfile_paired_flag_N+".tail");
	
	uint8_t chartoint[128];
	chartoint['0'] = 0;
	chartoint['1'] = 1;
	in_flag_N.close();
	in_flag_N.open(outfile_paired_flag_N);
	char chararray[8];
	uint8_t packedchar;
	for(uint64_t i = 0; i < numreads_clean/8; i++)
	{
		in_flag_N.read(chararray,8);
		
		packedchar = 128*chartoint[chararray[7]]+ 64*chartoint[chararray[6]]
					+ 32*chartoint[chararray[5]]+ 16*chartoint[chararray[4]]
					+  8*chartoint[chararray[3]]+  4*chartoint[chararray[2]]
					+  2*chartoint[chararray[1]]+  1*chartoint[chararray[0]];
		f_flag_N.write((char*)&packedchar,sizeof(uint8_t));
	}
	f_flag_N.close();
	in_flag_N.read(chararray,numreads_clean%8);
	for(int i=0; i<numreads_clean%8;i++)
		f_flag_N_tail << chararray[i];
	f_flag_N_tail.close();
	in_flag_N.close();
	remove((outfile_paired_flag_N).c_str());
	rename((outfile_paired_flag_N+".tmp").c_str(),(outfile_paired_flag_N).c_str());		

	return;
}
