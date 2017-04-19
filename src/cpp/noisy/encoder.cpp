#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>


long readlen;

std::string infile;
std::string infile_flag;
std::string infile_pos;
std::string infile_seq;
std::string outfile_seq;
std::string outfile_meta;
std::string outfile_pos;
std::string outfile_noise;
std::string outfile_noisepos;


char longtochar[] = {'A','C','G','T'};
long chartolong[128];
char enc_noise[128][128];

void encode();

std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos);

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos);

void getDataParams();

void setglobalarrays();

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/output/temp.dna";
	infile_pos = basedir + "/output/temppos.txt";
	infile_flag = basedir + "/output/tempflag.txt";
	
	outfile_seq = basedir + "/output/read_seq.txt";
	outfile_meta = basedir + "/output/read_meta.txt";
	outfile_pos = basedir + "/output/read_pos.txt";
	outfile_noise = basedir + "/output/read_noise.txt";
	outfile_noisepos = basedir + "/output/read_noisepos.txt";
	getDataParams(); //populate readlen
	setglobalarrays();
	encode();
	return 0;
}

void encode()
{
	std::cout<<"Encoding reads\n";
	std::ifstream f(infile);
	std::ifstream in_flag(infile_flag);
	std::ifstream in_pos(infile_pos);
	std::ofstream f_seq(outfile_seq);
	std::ofstream f_meta(outfile_meta);
	std::ofstream f_pos(outfile_pos);
	std::ofstream f_noise(outfile_noise);
	std::ofstream f_noisepos(outfile_noisepos);
	
	std::string current,ref;
	char c;
	std::vector<std::string> reads;
	std::vector<long> pos;
	long p;
	while(std::getline(f,current))
	{
		c = in_flag.get();
		in_pos>>p;
		if(c=='0'||reads.size()>10000000)//so that reads vector doesn't get too large
		{
			if(reads.size()!=0)
			{
				ref = buildcontig(reads,pos);
				writecontig(ref,pos,reads,f_seq,f_pos,f_noise,f_noisepos);
			}
			reads = {current};
			pos = {p};
		}
		else
		{
			reads.push_back(current);
			pos.push_back(p);
		}	
					
	}
	ref = buildcontig(reads,pos);
	writecontig(ref,pos,reads,f_seq,f_pos,f_noise,f_noisepos);
	f_meta << readlen << "\n";
	f.close();
	in_flag.close();
	in_pos.close();
	f_seq.close();
	f_meta.close();
	f_pos.close();
	f_noise.close();
	f_noisepos.close();
	std::cout << "Encoding done\n";
	return;
}


std::string buildcontig(std::vector<std::string> reads, std::vector<long> pos)
{
	if(reads.size() == 1)
		return reads[0];
	std::vector<std::array<long,4>> count(readlen,{0,0,0,0});
	for(long i = 0; i < readlen; i++)
		count[i][chartolong[reads[0][i]]] = 1;
	long prevpos = 0,currentpos;
	for(long j = 1; j < reads.size(); j++)
	{
		currentpos = prevpos + pos[j];
		count.insert(count.end(),pos[j],{0,0,0,0});
		for(long i = 0; i < readlen; i++)
			count[currentpos+i][chartolong[reads[j][i]]] += 1;
		prevpos = currentpos;
	}
	std::string ref(count.size(),'A');
	for(long i = 0; i < count.size(); i++)
	{
		long max = 0,indmax = 0;
		for(long j = 0; j < 4; j++)
			if(count[i][j]>max)
			{
				max = count[i][j];
				indmax = j;
			}
		ref[i] = longtochar[indmax];
	}
	return ref;
}

void writecontig(std::string ref,std::vector<long> pos, std::vector<std::string> reads, std::ofstream& f_seq, std::ofstream& f_pos, std::ofstream& f_noise, std::ofstream& f_noisepos)
{
	f_seq << ref;
	char c;
	if(reads.size() == 1)
	{
		f_noise << "\n";
		c = readlen;//(not pos[0] to handle breaks in read sequence due to limit on reads.size() - can't
				//assume  pos[0] = readlen)
		f_pos << c;
		return;
	}
	long prevj = 0;
	for(long j = 0; j < readlen; j++)
		if(reads[0][j] != ref[j])
		{
			f_noise<<enc_noise[ref[j]][reads[0][j]];
			c = j-prevj;
			f_noisepos<<c;
			prevj = j;
		}
	f_noise << "\n";
	c = readlen;// (to handle breaks in read sequence due to limit on reads.size()
	f_pos << c;
	long prevpos = 0,currentpos;
	for(long i = 1; i < reads.size(); i++)
	{
		currentpos = prevpos + pos[i];
		prevj = 0;
		for(long j = 0; j < readlen; j++)
			if(reads[i][j] != ref[currentpos+j])
			{
				f_noise<<enc_noise[ref[currentpos+j]][reads[i][j]];
				c = j-prevj;
				f_noisepos<<c;
				prevj = j;
			}
		f_noise << "\n";
		c = pos[i];
		f_pos << c;
		prevpos = currentpos;
	}
	return;
}

void setglobalarrays()
{
	//enc_noise uses substitution statistics from Minoche et al. 
	enc_noise['A']['C'] = '0';
	enc_noise['A']['G'] = '1';
	enc_noise['A']['T'] = '2';
	enc_noise['A']['N'] = '3';
	enc_noise['C']['A'] = '0';
	enc_noise['C']['G'] = '1';
	enc_noise['C']['T'] = '2';
	enc_noise['C']['N'] = '3';
	enc_noise['G']['T'] = '0';
	enc_noise['G']['A'] = '1';
	enc_noise['G']['C'] = '2';
	enc_noise['G']['N'] = '3';
	enc_noise['T']['G'] = '0';
	enc_noise['T']['C'] = '1';
	enc_noise['T']['A'] = '2';
	enc_noise['T']['N'] = '3';
	enc_noise['N']['A'] = '0';
	enc_noise['N']['G'] = '1';
	enc_noise['N']['C'] = '2';
	enc_noise['N']['T'] = '3';
	chartolong['A'] = 0;
	chartolong['C'] = 1;
	chartolong['G'] = 2;
	chartolong['T'] = 3;
	chartolong['N'] = 4;
	return;
}
	

void getDataParams()
{
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);

	std::getline(myfile, line);
	readlen = line.length();
	std::cout << "Read length: " << readlen << std::endl;
	myfile.close();
}
