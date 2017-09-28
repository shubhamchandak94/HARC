#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <bitset>
#include <cstdio>
#include <omp.h>
#include "config.h"

std::string outfile;
std::string outfile_clean;
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_order;
std::string infile_N;
std::string infile_order_N_pe;
std::string infile_singleton;
std::string infile_flag_N;


long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];
char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)


void decode();

void restore_order();

void restore_order_N();

void unpackbits();

std::bitset<2*readlen> stringtobitset(std::string s);

std::bitset<2*readlen> chartobitset(char *s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

void reverse_complement(char* s, char* s1);

std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);


void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

uint32_t numreads = 0;

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.tmp";
	outfile_clean = basedir + "/output/output_clean.dna";
	infile_seq = basedir + "/output/read_seq.txt";
	infile_meta = basedir + "/output/read_meta.txt";
	infile_pos = basedir + "/output/read_pos.txt";
	infile_noise = basedir + "/output/read_noise.txt";
	infile_noisepos = basedir + "/output/read_noisepos.txt";
	infile_rev = basedir + "/output/read_rev.txt";
	infile_order = basedir + "/output/read_order.bin";
	infile_N = basedir + "/output/input_N.dna";
	infile_order_N_pe =  basedir + "/output/read_order_N_pe.bin";
	infile_singleton = basedir + "/output/read_singleton.txt";
	infile_flag_N =  basedir + "/output/read_flag_N.txt";

//	getDataParams(); //populate readlen
	omp_set_num_threads(num_thr);
	setglobalarrays();
	decode();
	restore_order_N();
	restore_order();
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	unpackbits();
	uint32_t numreads_thr[num_thr];
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	numreads_thr[tid] = 0;
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f(outfile+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream f_seq(infile_seq+'.'+std::to_string(tid_e));
		std::ifstream f_pos(infile_pos+'.'+std::to_string(tid_e));
		std::ifstream f_noise(infile_noise+'.'+std::to_string(tid_e));
		std::ifstream f_noisepos(infile_noisepos+'.'+std::to_string(tid_e));
		std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e));
		std::ofstream f_N_tmp(infile_N+'.'+std::to_string(tid_e)+".tmp");
		std::ifstream f_flag_N(infile_flag_N+'.'+std::to_string(tid_e));

		char currentread[readlen+1],ref[readlen+1],revread[readlen+1];
		std::bitset<2*readlen> b;
		currentread[readlen] = '\0';
		revread[readlen] = '\0';
		ref[readlen] = '\0';
		std::string noise;
		char c;
		long pos;
		char flag_N; 

		while(f_pos >> std::noskipws >> c)//don't skip whitespaces
		{
			pos = (unsigned char)(c);
			if(pos!=0)
			{
				for(int i = 0; i <= readlen - 1 - pos;  i++)
					ref[i] = ref[i+pos];
				f_seq.get(ref+readlen-pos,pos+1);
			}
			strcpy(currentread,ref);	
			int prevnoisepos = 0,noisepos;
			std::getline(f_noise,noise);
			for(int i = 0; i < noise.size(); i++)
			{
				c = f_noisepos.get();
				noisepos = (unsigned char)c + prevnoisepos;
				currentread[noisepos] = dec_noise[ref[noisepos]][noise[i]];
				prevnoisepos = noisepos;	
			}
			c = f_rev.get();
			flag_N = f_flag_N.get();
			//if(strchr(currentread,'N')!=NULL)
			if(flag_N == '1')
			{
				if(c == 'd')
					f_N_tmp << currentread<<"\n";
				else
				{
					reverse_complement(currentread,revread);
					f_N_tmp << revread<<"\n";
				}
			}
			else
			{
				numreads_thr[tid]++;
				if(c == 'd')
				{
					b = chartobitset(currentread);
					f.write((char*)&b,sizeof(std::bitset<2*readlen>));
				}
				else
				{
					reverse_complement(currentread,revread);
					b = chartobitset(revread);
					f.write((char*)&b,sizeof(std::bitset<2*readlen>));
				}
			}
		}
		f.close();
		f_seq.close();
		f_pos.close();
		f_noise.close();
		f_noisepos.close();
		f_rev.close();
		f_N_tmp.close();
		f_flag_N.close();
	}
	}
	numreads = 0;
	for(int tid = 0; tid < num_thr; tid++)
		numreads += numreads_thr[tid];
	std::ofstream f(outfile,std::ios::binary);
	for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
	{
		std::ifstream f_in(outfile+'.'+std::to_string(tid_e),std::ios::binary);
		f << f_in.rdbuf();	
		f.clear();//clear error flags if rdbuf empty file	
		f_in.close();
	}
	std::ifstream f_singleton(infile_singleton);
	char currentread[readlen+1];
	std::bitset<2*readlen> b;
	currentread[readlen] = '\0';
	f_singleton.read(currentread,readlen);
	while(!f_singleton.eof())
	{	
		numreads++;
		b = chartobitset(currentread);
		f.write((char*)&b,sizeof(std::bitset<2*readlen>));
		f_singleton.read(currentread,readlen);
	}
	f_singleton.close();
	std::ofstream f_N_tmp(infile_N+".tmp",std::ios::binary);	
	for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
	{
		std::ifstream in_N_tmp(infile_N+'.'+std::to_string(tid_e)+".tmp");
		f_N_tmp << in_N_tmp.rdbuf();			
		f_N_tmp.clear();//clear error flags if rdbuf empty file	
		in_N_tmp.close();
	}
	std::ifstream f_N(infile_N);
	f_N_tmp << f_N.rdbuf();
	f_N_tmp.close();
	f_N.close();
	rename((infile_N+".tmp").c_str(), infile_N.c_str());
	std::cout<<"Decoding done\n";
	return;
}

void restore_order_N()
{
	std::string line;
	std::ifstream f_N(infile_N);
	uint32_t numreads_N = 0;
	while(std::getline(f_N,line))
		numreads_N++;
	f_N.close();
	f_N.open(infile_N);
		
	std::ifstream f_order(infile_order_N_pe,std::ios::binary);
	uint32_t *index_array = new uint32_t [numreads_N];
	char(*reads_N)[readlen+1] = new char [numreads_N][readlen+1];
	uint32_t order;
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		f_order.read((char*)&order,sizeof(uint32_t));
		index_array[order] = j;
		f_N.getline(reads_N[j],readlen+1);
	}
	f_N.close();
	std::ofstream out_N(infile_N);
	
	for(uint32_t j = 0; j < numreads_N; j++)
	{
		out_N << reads_N[index_array[j]] << "\n";
	}
	delete[] index_array;
	delete[] reads_N;
	f_order.close();
	out_N.close();
	return;
}

void restore_order()
{
	std::cout << "Restoring order\n";
	uint64_t max_bin_size;
	if(MAX_BIN_SIZE <= 3)
		max_bin_size = uint64_t(3)*200000000/7;
	else 
		max_bin_size = uint64_t(MAX_BIN_SIZE)*200000000/7;
	char s[readlen+1];
	s[readlen] = '\0';
	std::ofstream f_clean(outfile_clean);
	for (uint32_t i = 0; i <= numreads/max_bin_size; i++)
	{
		std::ifstream f_order(infile_order,std::ios::binary);
		std::ifstream f(outfile,std::ios::binary);
		auto numreads_bin = max_bin_size;
		if (i == numreads/max_bin_size)
			numreads_bin = numreads%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		std::bitset<2*readlen> *reads_bin = new std::bitset<2*readlen> [numreads_bin];	
		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f.seekg(uint64_t(j)*sizeof(std::bitset<2*readlen>), f.beg);
				f.read((char*)&reads_bin[pos],sizeof(std::bitset<2*readlen>));
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			bitsettostring(reads_bin[index_array[j]],s);
			f_clean << s << "\n";
		}
		delete[] index_array;
		delete[] reads_bin;
		f_order.close();
		f.close();
	}
	f_clean.close();
	return;
}

void unpackbits()
{
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f_seq(infile_seq+'.'+std::to_string(tid_e)+".tmp");
		std::ofstream f_rev(infile_rev+'.'+std::to_string(tid_e)+".tmp");
		std::ifstream in_seq(infile_seq+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream in_seq_tail(infile_seq+'.'+std::to_string(tid_e)+".tail");
		std::ifstream in_rev(infile_rev+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream in_rev_tail(infile_rev+'.'+std::to_string(tid_e)+".tail");
		char inttobase[4];
		inttobase[0] = 'A';
		inttobase[1] = 'C';
		inttobase[2] = 'G';
		inttobase[3] = 'T';
		
		uint8_t dnabin;
		in_seq.read((char*)&dnabin,sizeof(uint8_t));
		while(!in_seq.eof())
		{	
			f_seq << inttobase[dnabin%4];
			dnabin/=4; 	
			f_seq << inttobase[dnabin%4];
			dnabin/=4; 	
			f_seq << inttobase[dnabin%4];
			dnabin/=4; 	
			f_seq << inttobase[dnabin%4];
			in_seq.read((char*)&dnabin,sizeof(uint8_t));
		}
		in_seq.close();
		f_seq << in_seq_tail.rdbuf();
		in_seq_tail.close();
		
		//rev
		inttobase[0] = 'd';
		inttobase[1] = 'r';
		
		in_rev.read((char*)&dnabin,sizeof(uint8_t));
		while(!in_rev.eof())
		{	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			dnabin/=2; 	
			f_rev << inttobase[dnabin%2];
			in_rev.read((char*)&dnabin,sizeof(uint8_t));
		}
		in_rev.close();
		f_rev << in_rev_tail.rdbuf();
		in_rev_tail.close();
		
		f_seq.close();
		f_rev.close();
		remove((infile_seq+'.'+std::to_string(tid_e)).c_str());
		remove((infile_rev+'.'+std::to_string(tid_e)).c_str());
		remove((infile_seq+'.'+std::to_string(tid_e)+".tail").c_str());
		rename((infile_seq+'.'+std::to_string(tid_e)+".tmp").c_str(),(infile_seq+'.'+std::to_string(tid_e)).c_str());
		remove((infile_rev+'.'+std::to_string(tid_e)+".tail").c_str());
		rename((infile_rev+'.'+std::to_string(tid_e)+".tmp").c_str(),(infile_rev+'.'+std::to_string(tid_e)).c_str());	
	}//for end
	}//parallel end
	
	//singleton
	std::ifstream in_singleton(infile_singleton,std::ios::binary);
	std::ofstream f_singleton(infile_singleton+".tmp");
	std::ifstream in_singleton_tail(infile_singleton+".tail");
	char inttobase[4];
	inttobase[0] = 'A';
	inttobase[1] = 'C';
	inttobase[2] = 'G';
	inttobase[3] = 'T';
	
	uint8_t dnabin;
	in_singleton.read((char*)&dnabin,sizeof(uint8_t));
	while(!in_singleton.eof())
	{	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		dnabin/=4; 	
		f_singleton << inttobase[dnabin%4];
		in_singleton.read((char*)&dnabin,sizeof(uint8_t));
	}
	in_singleton.close();
	f_singleton << in_singleton_tail.rdbuf();
	in_singleton_tail.close();
	f_singleton.close();
	remove(infile_singleton.c_str());
	remove((infile_singleton+".tail").c_str());
	rename((infile_singleton+".tmp").c_str(),infile_singleton.c_str());		
	return;
}

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}


void setglobalarrays()
{
	dec_noise['A']['0'] = 'C';
	dec_noise['A']['1'] = 'G';
	dec_noise['A']['2'] = 'T';
	dec_noise['A']['3'] = 'N';
	dec_noise['C']['0'] = 'A';
	dec_noise['C']['1'] = 'G';
	dec_noise['C']['2'] = 'T';
	dec_noise['C']['3'] = 'N';
	dec_noise['G']['0'] = 'T';
	dec_noise['G']['1'] = 'A';
	dec_noise['G']['2'] = 'C';
	dec_noise['G']['3'] = 'N';
	dec_noise['T']['0'] = 'G';
	dec_noise['T']['1'] = 'C';
	dec_noise['T']['2'] = 'A';
	dec_noise['T']['3'] = 'N';
	dec_noise['N']['0'] = 'A';
	dec_noise['N']['1'] = 'G';
	dec_noise['N']['2'] = 'C';
	dec_noise['N']['3'] = 'T';
	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartorevchar['N'] = 'N';
	for(int i = 0; i < 64; i++)
		mask64[i] = 1;
	for(int i = 0; i < readlen; i++)
	{
		basemask[i]['A'][2*i] = 0;
		basemask[i]['A'][2*i+1] = 0;
		basemask[i]['C'][2*i] = 0;
		basemask[i]['C'][2*i+1] = 1;
		basemask[i]['G'][2*i] = 1;
		basemask[i]['G'][2*i+1] = 0;
		basemask[i]['T'][2*i] = 1;
		basemask[i]['T'][2*i+1] = 1;
		positionmask[i][2*i] = 1;
		positionmask[i][2*i+1] = 1;
	}		
	return;
}
	

std::bitset<2*readlen> stringtobitset(std::string s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	unsigned long long ull,rem;
	for(int i = 0; i < (2*readlen)/64+1; i++)
	{	
		ull = (b&mask64).to_ullong();
		b>>=64;
		for(int j = 32*i  ; j < 32*i+32 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%4];	
			ull/=4;
		}
	}
	return;
}

std::bitset<2*readlen> chartobitset(char *s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

