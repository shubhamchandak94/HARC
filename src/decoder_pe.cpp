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

std::string outfile_1;
std::string outfile_2;
std::string outfile_order_2;
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_N;
std::string infile_singleton;
std::string infilenumreads;
std::string infile_order_paired;
std::string infile_paired_flag_first;


int readlen, MAX_BIN_SIZE, num_thr, num_thr_e;
uint32_t numreads, numreads_by_2; 

typedef std::bitset<3*MAX_READ_LEN> bitset;

long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];
char revinttochar[8] = {'A','N','G','#','C','#','T','#'};//used in bitsettostring
bitset basemask[MAX_READ_LEN][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
bitset positionmask[MAX_READ_LEN];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
bitset mask63;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)


void decode(bool *flag_first);

void generate_order_from_paired(bool *flag_first);

void restore_order();

void unpackbits();

bitset stringtobitset(std::string s);

bitset chartobitset(char *s);

void bitsettostring(bitset b,char *s);

void reverse_complement(char* s, char* s1);

std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile_1 = basedir + "/output_1.dna";
	outfile_2 = basedir + "/output_2.dna";
	outfile_order_2 = basedir + "/read_order_2.bin";
	infile_seq = basedir + "/read_seq.txt";
	infile_meta = basedir + "/read_meta.txt";
	infile_pos = basedir + "/read_pos.txt";
	infile_noise = basedir + "/read_noise.txt";
	infile_noisepos = basedir + "/read_noisepos.txt";
	infile_rev = basedir + "/read_rev.txt";
	infile_N = basedir + "/unaligned_N.txt";
	infile_singleton = basedir + "/unaligned_singleton.txt";
	infilenumreads = basedir + "/numreads.bin";
	infile_order_paired = basedir + "/read_order_paired.bin";
	infile_paired_flag_first = basedir + "/read_paired_flag_first.bin";
	
	readlen = atoi(argv[2]);
	MAX_BIN_SIZE = atoi(argv[3]);
	num_thr = atoi(argv[4]);
	num_thr_e = atoi(argv[5]);

	std::ifstream f_numreads(infilenumreads, std::ios::binary);
	f_numreads.seekg(4);
	f_numreads.read((char*)&numreads,sizeof(uint32_t));
	f_numreads.close();
	numreads_by_2 = numreads/2;

	omp_set_num_threads(num_thr);
	setglobalarrays();
	
	bool *flag_first = new bool [numreads];//flag indicating first reads of pair
	
	generate_order_from_paired(flag_first);
	
	decode(flag_first);
	restore_order();
	return 0;
}

void decode(bool *flag_first)
{
	std::cout << "Decoding reads\n";
	unpackbits();
	uint32_t numreads_thr[num_thr];
	//first calculate numreads in each thread, this is needed for accessing flag_first to know which reads are first in pair
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		numreads_thr[tid] = 0;
		for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
		{
			std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e), std::ifstream::ate | std::ifstream::binary);	
			numreads_thr[tid] += f_rev.tellg();//size of f_rev file
			f_rev.close();
		}
	}		
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t pos_in_flag_first = 0;
	if(tid != 0)
	{
		for(int i = 0; i < tid; i++)
			pos_in_flag_first += numreads_thr[i];
	}
	
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f_1(outfile_1+'.'+std::to_string(tid_e));
		std::ofstream f_2(outfile_2+'.'+std::to_string(tid_e),std::ios::binary);
		std::ifstream f_seq(infile_seq+'.'+std::to_string(tid_e));
		std::ifstream f_pos(infile_pos+'.'+std::to_string(tid_e));
		std::ifstream f_noise(infile_noise+'.'+std::to_string(tid_e));
		std::ifstream f_noisepos(infile_noisepos+'.'+std::to_string(tid_e));
		std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e));
	
		char currentread[MAX_READ_LEN+1],ref[MAX_READ_LEN+1],revread[MAX_READ_LEN+1];
		bitset b;
		currentread[readlen] = '\0';
		revread[readlen] = '\0';
		ref[readlen] = '\0';
		std::string noise;
		char c;
		long pos; 

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
			if(c == 'd')
			{
				if(flag_first[pos_in_flag_first])
				{
					f_1 << currentread << "\n";
				}	
				else
				{
					b = chartobitset(currentread);
					f_2.write((char*)&b,sizeof(bitset));
				}	
			}
			else
			{
				reverse_complement(currentread,revread);
				if(flag_first[pos_in_flag_first])
				{
					f_1 << revread << "\n";
				}	
				else
				{
					b = chartobitset(revread);
					f_2.write((char*)&b,sizeof(bitset));
				}	
			}
			pos_in_flag_first++;
		}
		f_1.close();
		f_2.close();
		f_seq.close();
		f_pos.close();
		f_noise.close();
		f_noisepos.close();
		f_rev.close();
	}
	}
	std::ofstream f_1(outfile_1);
	std::ofstream f_2(outfile_2,std::ios::binary);
	for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
	{
		std::ifstream f_in_1(outfile_1+'.'+std::to_string(tid_e));
		std::ifstream f_in_2(outfile_2+'.'+std::to_string(tid_e),std::ios::binary);
		f_1 << f_in_1.rdbuf();
		f_2 << f_in_2.rdbuf();
		f_1.clear();//clear error flags if f_in is empty	
		f_2.clear();//clear error flags if f_in is empty	
		f_in_1.close();
		f_in_2.close();
	}
	uint32_t pos_in_flag_first = 0;
	for(int i = 0; i < num_thr; i++)
		pos_in_flag_first += numreads_thr[i];
	std::ifstream f_singleton(infile_singleton);
	char currentread[MAX_READ_LEN+1];
	bitset b;
	f_singleton.read(currentread,readlen);
	while(!f_singleton.eof())
	{	
		if(flag_first[pos_in_flag_first])
		{
			f_1 << currentread << "\n";
		}
		else
		{			
			b = chartobitset(currentread);
			f_2.write((char*)&b,sizeof(bitset));
		}
		f_singleton.read(currentread,readlen);
		pos_in_flag_first++;
	}
	f_singleton.close();
	std::ifstream f_N(infile_N);
	f_N.read(currentread,readlen);
	while(!f_N.eof())
	{
		if(flag_first[pos_in_flag_first])
		{
			f_1 << currentread << "\n";
		}
		else
		{			
			b = chartobitset(currentread);
			f_2.write((char*)&b,sizeof(bitset));
		}
		f_N.read(currentread,readlen);
		pos_in_flag_first++;
	}
	f_N.close();
	f_1.close();
	f_2.close();
	std::cout<<"Decoding done\n";
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
	char s[MAX_READ_LEN+1];
	std::ofstream fout(outfile_2+".tmp");
	for (uint32_t i = 0; i <= numreads_by_2/max_bin_size; i++)
	{
		std::ifstream f_order(outfile_order_2,std::ios::binary);
		std::ifstream f(outfile_2,std::ios::binary);
		auto numreads_bin = max_bin_size;
		if (i == numreads_by_2/max_bin_size)
			numreads_bin = numreads_by_2%max_bin_size;
		uint32_t *index_array = new uint32_t [numreads_bin];
		bitset *reads_bin = new bitset [numreads_bin];	
		uint32_t order,pos = 0;
		for(uint32_t j = 0; j < numreads_by_2; j++)
		{
			f_order.read((char*)&order,sizeof(uint32_t));
			if (order >= i*max_bin_size && order < i*max_bin_size + numreads_bin)
			{
				index_array[order-i*max_bin_size] = pos;
				f.seekg(uint64_t(j)*sizeof(bitset), f.beg);
				f.read((char*)&reads_bin[pos],sizeof(bitset));
				pos++;
			}
		}
		for(uint32_t j = 0; j < numreads_bin; j++)
		{
			bitsettostring(reads_bin[index_array[j]],s);
			fout << s << "\n";
		}
		delete[] index_array;
		delete[] reads_bin;
		f_order.close();
		f.close();
	}
	fout.close();
	remove(outfile_2.c_str());
	rename((outfile_2+".tmp").c_str(),outfile_2.c_str());
	return;
}

void generate_order_from_paired(bool *flag_first)
{
	uint32_t *read_order = new uint32_t[numreads];
	std::fill(read_order, read_order+numreads, numreads);
	std::ifstream in_order_paired(infile_order_paired,std::ios::binary);
	std::ifstream in_paired_flag_first(infile_paired_flag_first);
	
	char c;
	uint32_t order_paired;
	uint32_t current_pair = 0;
	//decode order_paired and flag_first into read_order 
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(read_order[i] != numreads)//this position already filled
			continue;
		in_paired_flag_first.get(c);
		in_order_paired.read((char*)&order_paired, sizeof(uint32_t));
		if(c == '1')
		{
			read_order[i] = current_pair;
			flag_first[i] = 1;
			read_order[i+order_paired] = current_pair+numreads_by_2;
			flag_first[i+order_paired] = 0;
		}
		else
		{
			read_order[i] = current_pair+numreads_by_2;
			flag_first[i] = 0;
			read_order[i+order_paired] = current_pair;
			flag_first[i+order_paired] = 1;
		}
		current_pair++;	
	}
	
	//now we'll need the inverse index, before doing that write to file to save memory
	std::ofstream f_temp(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		f_temp.write((char*)&read_order[i], sizeof(uint32_t));
	}
	f_temp.close();
	
	//build inverse index
	std::ifstream fin_temp(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		read_order[order_paired] = i;
	}
	fin_temp.close();
	
	//now, go through the reads in order, and for each read which is 1st in paired end, 
	//store the location of its pair 
	f_temp.open(outfile_order_2+".tmp",std::ios::binary);
	fin_temp.open(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		if(flag_first[i])
		{
			uint32_t pos_of_pair = read_order[order_paired+numreads_by_2];
			f_temp.write((char*)&pos_of_pair, sizeof(uint32_t));
		}
	}
	fin_temp.close();
	f_temp.close();
	
	//store the inverse index for the 2nd reads of pair
	fin_temp.open(outfile_order_2+".tmp",std::ios::binary);
	for(uint32_t i = 0; i < numreads_by_2; i++)
	{
		fin_temp.read((char*)&order_paired, sizeof(uint32_t));
		read_order[order_paired] = i;
	}
	fin_temp.close();
	remove((outfile_order_2+".tmp").c_str());
	
	//Finally, write the order for second reads of pair 
	f_temp.open(outfile_order_2,std::ios::binary);
	for(uint32_t i = 0; i < numreads; i++)
	{
		if(!flag_first[i])
		{
			f_temp.write((char*)&read_order[i], sizeof(uint32_t));
		}
	}
	f_temp.close();
	delete[] read_order;
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
	
	//flag_first
	std::ifstream in_flag_first(infile_paired_flag_first,std::ios::binary);
	std::ofstream f_flag_first(infile_paired_flag_first+".tmp");
	std::ifstream in_flag_first_tail(infile_paired_flag_first+".tail");
	char inttochar[2];
	inttochar[0] = '0';
	inttochar[1] = '1';
	
	uint8_t packedchar;
	in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	while(!in_flag_first.eof())
	{	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		f_flag_first << inttochar[packedchar%2];
		packedchar/=2; 	
		in_flag_first.read((char*)&packedchar,sizeof(uint8_t));
	}
	in_flag_first.close();
	f_flag_first << in_flag_first_tail.rdbuf();
	in_flag_first_tail.close();
	f_flag_first.close();
	remove(infile_paired_flag_first.c_str());
	remove((infile_paired_flag_first+".tail").c_str());
	rename((infile_paired_flag_first+".tmp").c_str(),infile_paired_flag_first.c_str());
	return;
}

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	s1[readlen] = '\0';
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
	for(int i = 0; i < 63; i++)
		mask63[i] = 1;
	for(int i = 0; i < readlen; i++)
	{
		basemask[i]['A'][3*i] = 0;
		basemask[i]['A'][3*i+1] = 0;
		basemask[i]['A'][3*i+2] = 0;
		basemask[i]['C'][3*i] = 0;
		basemask[i]['C'][3*i+1] = 0;
		basemask[i]['C'][3*i+2] = 1;
		basemask[i]['G'][3*i] = 0;
		basemask[i]['G'][3*i+1] = 1;
		basemask[i]['G'][3*i+2] = 0;
		basemask[i]['T'][3*i] = 0;
		basemask[i]['T'][3*i+1] = 1;
		basemask[i]['T'][3*i+2] = 1;
		basemask[i]['N'][3*i] = 1;
		basemask[i]['N'][3*i+1] = 0;
		basemask[i]['N'][3*i+2] = 0;
		positionmask[i][3*i] = 1;
		positionmask[i][3*i+1] = 1;
		positionmask[i][3*i+2] = 1;
	}		
	return;
}
	

bitset stringtobitset(std::string s)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void bitsettostring(bitset b,char *s)
{
	unsigned long long ull,rem;
	for(int i = 0; i < (3*readlen)/63+1; i++)
	{	
		ull = (b&mask63).to_ullong();
		b>>=63;
		for(int j = 21*i  ; j < 21*i+21 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%8];	
			ull/=8;
		}
	}
	s[readlen] = '\0';
	return;
}

bitset chartobitset(char *s)
{
	bitset b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}
