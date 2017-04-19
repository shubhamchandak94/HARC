#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>


long readlen;

std::string outfile;
std::string infile_seq;
std::string infile_meta;
std::string infile_pos;
std::string infile_noise;
std::string infile_noisepos;
std::string infile_rev;
std::string infile_N;


char longtochar[] = {'A','C','G','T'};
long chartolong[128];
char dec_noise[128][128];
char chartorevchar[128];

void decode();

void reverse_complement(char* s, char* s1);
std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	outfile = basedir + "/output.dna";
	
	infile_seq = basedir + "/output/read_seq.txt";
	infile_meta = basedir + "/output/read_meta.txt";
	infile_pos = basedir + "/output/read_pos.txt";
	infile_noise = basedir + "/output/read_noise.txt";
	infile_noisepos = basedir + "/output/read_noisepos.txt";
	infile_rev = basedir + "/output/read_rev.txt";
	infile_N = basedir + "/output/input_N.dna";
	getDataParams(); //populate readlen
	setglobalarrays();
	decode();
	return 0;
}

void decode()
{
	std::cout << "Decoding reads\n";
	std::ofstream f(outfile);
	std::ifstream f_seq(infile_seq);
	std::ifstream f_pos(infile_pos);
	std::ifstream f_noise(infile_noise);
	std::ifstream f_noisepos(infile_noisepos);
	std::ifstream f_rev(infile_rev);
	std::ifstream f_N(infile_N);
	
	char currentread[readlen+1],ref[readlen+1],revread[readlen+1];
	currentread[readlen] = '\0';
	revread[readlen] = '\0';
	ref[readlen] = '\0';
	std::string noise;
	char c;
	long pos; 
	while(f_pos >> std::noskipws >> c)//don't skip whitespaces
	{
		pos = c;
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
			noisepos = c + prevnoisepos;
			currentread[noisepos] = dec_noise[ref[noisepos]][noise[i]];
			prevnoisepos = noisepos;	
		}
		c = f_rev.get();
		if(c == 'd')
			f << currentread<<"\n";
		else
		{
			reverse_complement(currentread,revread);
			f << revread<<"\n";
		}
	}
	f << f_N.rdbuf();		
	f.close();
	f_seq.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	f_rev.close();
	f_N.close();
	std::cout<<"Decoding done\n";
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


