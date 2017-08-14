#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

long readlen;
int num_thr, num_thr_e; //number of decoding and encoding threads

std::string outfile;
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_N;
std::string infile_singleton;


long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];

void decode();

void unpackbits();

void reverse_complement(char* s, char* s1);
std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output/output.dna";
	num_thr = atoi(argv[2]);
	num_thr_e = atoi(argv[3]);
	omp_set_num_threads(num_thr);
	
	infile_seq = basedir + "/output/read_seq.txt";
	infile_meta = basedir + "/output/read_meta.txt";
	infile_pos = basedir + "/output/read_pos.txt";
	infile_noise = basedir + "/output/read_noise.txt";
	infile_noisepos = basedir + "/output/read_noisepos.txt";
	infile_rev = basedir + "/output/read_rev.txt";
	infile_N = basedir + "/output/input_N.dna";
	infile_singleton = basedir + "/output/read_singleton.txt";
	getDataParams(); //populate readlen
	setglobalarrays();
	decode();
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	unpackbits();

	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	for(int tid_e = tid*num_thr_e/num_thr; tid_e < (tid+1)*num_thr_e/num_thr; tid_e++)
	{
		std::ofstream f(outfile+'.'+std::to_string(tid_e));
		std::ifstream f_seq(infile_seq+'.'+std::to_string(tid_e));
		std::ifstream f_pos(infile_pos+'.'+std::to_string(tid_e));
		std::ifstream f_noise(infile_noise+'.'+std::to_string(tid_e));
		std::ifstream f_noisepos(infile_noisepos+'.'+std::to_string(tid_e));
		std::ifstream f_rev(infile_rev+'.'+std::to_string(tid_e));
		std::ofstream f_N_tmp(infile_N+'.'+std::to_string(tid_e)+".tmp");
		
		char currentread[readlen+1],ref[readlen+1],revread[readlen+1];
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
				noisepos = (unsigned char)(c) + prevnoisepos;
				currentread[noisepos] = dec_noise[ref[noisepos]][noise[i]];
				prevnoisepos = noisepos;	
			}
			c = f_rev.get();
			if(strchr(currentread,'N')!=NULL)
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
				if(c == 'd')
					f << currentread<<"\n";
				else
				{
					reverse_complement(currentread,revread);
					f << revread<<"\n";
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
	}//for end
	}//parallel end
	std::ofstream f(outfile);
	for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
	{
		std::ifstream f_in(outfile+'.'+std::to_string(tid_e));
		f << f_in.rdbuf();			
		f_in.close();
	}
	std::ifstream f_singleton(infile_singleton);
	char currentread[readlen+1];
	currentread[readlen] = '\0';
	f_singleton.read(currentread,readlen);
	while(!f_singleton.eof())
	{	
		f << currentread << "\n";
		f_singleton.read(currentread,readlen);	
	}
	f_singleton.close();	
	for(int tid_e = 0; tid_e < num_thr_e; tid_e++)
	{
		std::ifstream f_N(infile_N+'.'+std::to_string(tid_e)+".tmp");
		f << f_N.rdbuf();			
		f_N.close();
	}
	std::ifstream f_N(infile_N);
	f << f_N.rdbuf();
	f_N.close();
	std::cout<<"Decoding done\n";
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
	return;
}
	

void getDataParams()
{
	std::string line;
	std::ifstream myfile(infile_meta, std::ifstream::in);
	
	myfile >> readlen;
	std::cout << "Read length: " << readlen << std::endl;
	myfile.close();
}
