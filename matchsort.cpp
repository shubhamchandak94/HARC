#include <iostream>
#include <fstream>
#include <bitset>
#include <unordered_set>
#include <sparsehash/sparse_hash_map>
#include <string>
#include <vector>

#define infile "chrom22_50x_noRC_random.dna"
#define outfile "temp.dna"
#define readlen 100
#define matchlen 20//has significant effect on size of dictionary
#define maxmatch 20
#define numreads 17500000


std::vector<bool> stringtobitset(std::string s)
{
	std::vector<bool> b(2*readlen,0);
	for(int i = 0; i < readlen; i++)//reverse so bitset is in correct order
	{	
		switch(s[i])
		{
			case 'A':	b[2*i] = 0;
					b[2*i+1] = 0;
					break;
			case 'C':	b[2*i] = 0;
					b[2*i+1] = 1;
					break;
			case 'G':	b[2*i] = 1;
					b[2*i+1] = 0;
					break;
			case 'T':	b[2*i] = 1;
					b[2*i+1] = 1;
					break;
		}
	}
	return b;
}
					
					

void readDnaFile(std::vector<bool> *read);

void constructdictionary(std::bitset<2*readlen> *read, google::sparse_hash_map<std::string,std::vector<int>> &dict);

void reorder();

void writetofile();


int main()
{
	std::vector<bool> *read = new std::vector<bool> [numreads];
	std::cout << "Reading file\n";
	readDnaFile(read);
	google::sparse_hash_map<std::string,std::vector<int>> dict;
	//using vector instead of list to save some space (at the cost of linear time required to delete elements)
	std::cout << "Constructing dictionary\n";
//	constructdictionary(read,dict);
	return 0;
}


void readDnaFile(std::vector<bool> *read)
{
	std::ifstream f(infile, std::ifstream::in);
	std::string s;
	for(int i = 0; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads read"<<"\n";
		f >> s;
		read[i] = stringtobitset(s);
	}
	f.close();
	return;
}
		
/*
void constructdictionary(std::bitset<2*readlen> *read, google::sparse_hash_map<std::string,std::vector<int>> &dict)
{
	for(int i = 0; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout << i/1000000 << "M reads read"<<"\n";
		std::string s = read[i].to_string().substr(0,2*matchlen);
		if(dict.count(s) == 1)
			dict[s].push_back(i);
		else
		 	dict[s] = {i};
	}

	return;
}
*/
