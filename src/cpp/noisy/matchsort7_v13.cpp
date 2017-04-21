#include <iostream>
#include <fstream>
#include <bitset>
#include <sparsepp/spp.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include "config.h"
#include <omp.h>

//#define readlen 100


uint32_t numreads = 0;

std::string infile;
std::string outfile;
std::string outfileRC;
std::string outfileflag;
std::string outfilepos;


void generateindexmasks(std::bitset<2*readlen> *mask1)
//masks for dictionary positions
{
	for(int i = 0; i < numdict; i++)
		mask1[i].reset();
	for(int i = 2*dict1_start; i < 2*(dict1_end+1); i++)
		mask1[0][i] = 1;
	for(int i = 2*dict2_start; i < 2*(dict2_end+1); i++)
		mask1[1][i] = 1;
	return;
}

char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];
int chartoint[128];
char chartobit[128][2];
int charinttoint[128];
std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)

std::bitset<2*readlen> stringtobitset(std::string s);

void bitsettostring(std::bitset<2*readlen> b,char *s);

void readDnaFile(std::bitset<2*readlen> *read);

void getDataParams();

void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict);

void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask);
//mask for zeroing the end bits (needed while reordering to compute Hamming distance between shifted reads)

uint32_t findread(std::vector<uint32_t> &s, std::bitset<2*readlen> *read, std::bitset<2*readlen> &b, std::bitset<2*readlen> &mask);
//search for read matching current read

uint32_t findread_parallel(std::vector<uint32_t> &s, std::bitset<2*readlen> *read, std::bitset<2*readlen> &b, std::bitset<2*readlen> &mask);

std::vector<uint32_t> dict_union(std::bitset<2*readlen> &b1, spp::sparse_hash_map<uint64_t,uint32_t*> *dict, int* dict_start, std::bitset<2*readlen>* mask);
//take union of reads in different dictionary bins

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict,std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool>& flagvec, std::vector<uint32_t>& readpos);

void writetofile(std::bitset<2*readlen> *read, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos);

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift);
//update the clean reference and count matrix

void reverse_complement(char *s, char *s1);

std::bitset<2*readlen> chartobitset(char &s);

void setglobalarrays();
//setting arrays like chartoint etc.

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/input_clean.dna";
	outfile = basedir + "/output/temp.dna";
	outfileRC = basedir + "/output/read_rev.txt";
	outfileflag = basedir + "/output/tempflag.txt";
	outfilepos = basedir + "/output/temppos.txt";
	getDataParams(); //populate numreads, readlen
	
	setglobalarrays();
	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	std::cout << "Reading file: " << infile << std::endl;
	readDnaFile(read);
	std::cout << "Constructing dictionaries\n";
	spp::sparse_hash_map<uint64_t,uint32_t*> dict[numdict];
	constructdictionary(read,dict);
	std::vector<uint32_t> sortedorder(numreads),readpos(numreads);
	std::vector<bool> revcomp(numreads),flagvec(numreads);
	std::cout << "Reordering reads\n";
	reorder(read,dict,sortedorder,revcomp,flagvec,readpos);
	std::cout << "Writing to file\n";
	writetofile(read,sortedorder,revcomp,flagvec,readpos);	
	delete[] read;
	std::cout << "Done!\n";
	return 0;
}

void setglobalarrays()
{
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartoint['A'] = 0;
	chartoint['C'] = 1;
	chartoint['G'] = 2;
	chartoint['T'] = 3;
	charinttoint['0'] = 0;
	charinttoint['1'] = 1;
	chartobit['A'][0] = '0';
	chartobit['A'][1] = '0';
	chartobit['C'][0] = '0';
	chartobit['C'][1] = '1';
	chartobit['G'][0] = '1';
	chartobit['G'][1] = '0';
	chartobit['T'][0] = '1';
	chartobit['T'][1] = '1';
	for(int i = 0; i < readlen; i++)
	{
		if(i < 64)
			mask64[i] = 1;
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
	
/*
std::bitset<2*readlen> stringtobitset(std::string s)
{
	char *s2 = &(s[0]);
	char s1[2*readlen+1];
	for(int i = 0; i < readlen; i++)
	{
		s1[2*(readlen-i-1)+1] = chartobit[s2[i]][0];
		s1[2*(readlen-i-1)] = chartobit[s2[i]][1];
	}
	s1[2*readlen] = '\0';
	std::bitset<2*readlen> b(s1);
	return b;
}
*/

std::bitset<2*readlen> stringtobitset(std::string s)
{
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void getDataParams()
{
	uint32_t number_of_lines = 0;
	std::string line;
	std::ifstream myfile(infile, std::ifstream::in);

	std::getline(myfile, line);
	int read_length = line.length();
	myfile.seekg(0, myfile.beg);

	while (std::getline(myfile, line))
	{
	++number_of_lines;
	int _len = line.length(); 
	if( _len != read_length) 
	{
		std::cout << "Read lengths are not the same!: " << read_length << " , " << _len << std::endl;
	}
	}
	//readlen = read_length;
	numreads = number_of_lines;
	std::cout << "Read length: " << read_length << std::endl;
	std::cout << "Number of reads: " << number_of_lines << std::endl;
	myfile.close();
}


void readDnaFile(std::bitset<2*readlen> *read)
{
	std::ifstream f(infile, std::ifstream::in);
	f.seekg(0, f.beg);
	std::string s;
	for(uint32_t i = 0; i < numreads; i++)
	{
		f >> s;
		read[i] = stringtobitset(s);
	}
	f.close();
	return;
}
		

void constructdictionary(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict)
{
	std::bitset<2*readlen> b;
	uint64_t ull;
	std::bitset<2*readlen> mask[numdict];
	generateindexmasks(mask);
	int dict_start[2] = {dict1_start,dict2_start};
	for(int j = 0; j < numdict; j++)
	{
		//find number of times each key occurs
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&mask[j];
			ull = (b>>2*dict_start[j]).to_ullong();
			if(dict[j].count(ull) == 1)
				(*dict[j][ull])++;
			else
			{
				dict[j][ull]=new uint32_t;
				(*dict[j][ull]) = 1;
			}
		}
		//allocate memory for each bin (number of reads with given key + 1) 1 for storing the length
		for(auto it = dict[j].begin(); it !=  dict[j].end(); ++it)
		{
			uint32_t binsize = *(it->second);
			delete it->second;
			dict[j][it->first] =  new uint32_t[binsize+1];
			dict[j][it->first][0] = 1;
		}
		//fill in the read ids in each bin, dict[j][ull][0] stores the position where next id is put - at the
		//end it stores the size of the array
		for(uint32_t i = 0; i < numreads; i++)
		{
			b = read[i]&mask[j];
			ull = (b>>2*dict_start[j]).to_ullong();
			dict[j][ull][dict[j][ull][0]++] = i;
		}

	}
	return;
}
/*
uint32_t findread_parallel(std::vector<uint32_t> &s, std::bitset<2*readlen> *read, std::bitset<2*readlen> &b, std::bitset<2*readlen> &mask)
//based on http://stackoverflow.com/questions/9793791/parallel-openmp-loop-with-break-statement
{
	auto N = s.size();
	if(N < 100)
		return findread(s,read,b,mask);
	bool done = false;
	uint32_t give = 0,k = numreads;
	#pragma omp parallel
	{
		uint32_t i, stop;
		#pragma omp critical
		{
			i = give;
			give += N/omp_get_num_threads();
			stop = give;
		
	
			if(omp_get_thread_num() == omp_get_num_threads()-1)
				stop = N;
		}
		while(i<stop && !done)
		{
			if((b^(read[s[i]]&mask)).count()<=thresh)
			{
				#pragma omp critical
				if(done == false)
				{
					done = true;
					k = s[i];
				}			
			}
			i++;
		}
	}
	return k;
}
*/

uint32_t findread(std::vector<uint32_t> &s, std::bitset<2*readlen> *read, std::bitset<2*readlen> &b, std::bitset<2*readlen> &mask)
{
	for (auto it = s.rbegin() ; it != s.rend(); ++it)
		if((b^(read[*it]&mask)).count()<=thresh)
			return *it;
	return numreads;
}

uint32_t findread_array(uint32_t* s, std::bitset<2*readlen> *read, std::bitset<2*readlen> &b, std::bitset<2*readlen> &mask)
{
	if(s[0] <= 100)
	{
		for (uint32_t i = s[0]-1 ; i >= 1; i--)
			if((b^(read[s[i]]&mask)).count()<=thresh)
				return s[i];
		return numreads;
	}
	else
	{
		auto N = s[0]-1;
		bool done = false;
		uint32_t give = 1,k = numreads;
		#pragma omp parallel
		{
			uint32_t i, stop;
			#pragma omp critical
			{
				i = give;
				give += N/omp_get_num_threads();
				stop = give;
			
		
				if(omp_get_thread_num() == omp_get_num_threads()-1)
					stop = N;
			}
			while(i<stop && !done)
			{
				if((b^(read[s[i]]&mask)).count()<=thresh)
				{
					#pragma omp critical
					if(done == false)
					{
						done = true;
						k = s[i];
					}			
				}
				i++;
			}
		}
		return k;
	}
}

std::vector<uint32_t> dict_union(std::bitset<2*readlen> &b1, spp::sparse_hash_map<uint64_t,uint32_t*> *dict, int* dict_start, std::bitset<2*readlen>* mask)
{
	uint64_t ull;
	std::vector<uint32_t> s;
	std::bitset<2*readlen> b;
	//take union of reads in different dictionary bins
	for(int l = 0; l < numdict; l++)
	{
		b = b1&mask[l];
		ull = (b>>2*dict_start[l]).to_ullong();
		if(dict[l].count(ull) == 1)
		{
			auto dictbin = dict[l][ull];
			if(s.empty())
				s.insert(s.end(),dictbin+1,dictbin+dictbin[0]);
			else
			{
				std::vector<uint32_t> temp;
				std::set_union(s.begin(),s.end(),dictbin+1,dictbin+dictbin[0],std::back_inserter(temp));
				s = temp;
			}
		}
	}
	return s;
}

void reorder(std::bitset<2*readlen> *read, spp::sparse_hash_map<uint64_t,uint32_t*> *dict, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos)
{	
	std::bitset<2*readlen> mask[maxmatch];
	std::bitset<2*readlen> revmask[maxmatch];
	std::bitset<2*readlen> ref,revref;
	int count[4][readlen];
	generatemasks(mask,revmask);
	std::bitset<2*readlen> mask1[numdict];
	generateindexmasks(mask1);
	bool *remainingreads = new bool[numreads];
	std::fill(remainingreads, remainingreads+numreads,1);
	uint32_t remainingpos = numreads-2;//used for searching next unmatched read when no match is found
	//we go through remainingreads array from behind as that speeds up deletion from bin arrays
	int dict_start[2] = {dict1_start,dict2_start};
	uint32_t unmatched = 0;
	uint32_t current = numreads-1;//picking last read as that speeds up deletion from bin arrays
	bool flag = 0;
	//flag to check if match was found or not
	sortedorder[0] = current;
	revcomp[0] = 0;//for direct
	flagvec[0] = 0;//for unmatched
	readpos[0] = readlen;//shift wrt previous read (always readlen for unmatched reads)
	updaterefcount(read[current],ref,revref,count,true,false,0);
	std::bitset<2*readlen> b;
	uint64_t ull;
	for(uint32_t i = 1; i < numreads; i++)
	{
		if(i%1000000 == 0)
			std::cout<< i/1000000 << "M reads done. \n";
		remainingreads[current] = 0;
		//delete the read from the corresponding dictionary bins
		for(int l = 0; l < numdict; l++)
		{	b = read[current]&mask1[l];
			ull = (b>>2*dict_start[l]).to_ullong();
			auto dictbin = dict[l][ull];
			if(flag==1)//if unmatched, we always get last read in bin (so no need to shift)
			{
				uint32_t pos = std::lower_bound(dictbin+1,dictbin+dictbin[0],current)-dictbin;
			//binary search since dict[l][ull] is sorted array
				std::move(dictbin+pos+1,dictbin+dictbin[0],dictbin+pos);
			}
			dictbin[0]--;//decrement number of elements
			if(dictbin[0] == 1)//empty (1 to store the size)
			{
				delete[] dictbin;
				dict[l].erase(ull);
			}
		}
			
		
		flag = 0;
		for(int j = 0; j < maxmatch; j++)
		{
/*
			std::vector<uint32_t> s = dict_union(ref,dict,dict_start,mask1);
			if(!s.empty())
			{
				uint32_t k = findread(s,read,ref,mask[j]);
				if(k != numreads)
				{
					current = k;
					flag = 1;
					updaterefcount(read[current],ref,revref,count,false,false,j);
					revcomp[i] = 0;
					sortedorder[i] = current;
					flagvec[i] = 1;//for matched
					readpos[i] = j;
					break;
				}	
			}
*/		
			for(int l = 0; l < numdict; l++)
			{
				b = ref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				if(dict[l].count(ull) == 1)
				{
					uint32_t k = findread_array(dict[l][ull],read,ref,mask[j]);
					if(k != numreads)
					{
						current = k;
						flag = 1;
						updaterefcount(read[current],ref,revref,count,false,false,j);
						revcomp[i] = 0;
						sortedorder[i] = current;
						flagvec[i] = 1;//for matched
						readpos[i] = j;
						break;
					}
				}
			}
			if(flag==1)
				break;	

			ref>>=2;

			for(int l = 0; l < numdict; l++)
			{
				b = revref&mask1[l];
				ull = (b>>2*dict_start[l]).to_ullong();
				if(dict[l].count(ull) == 1)
				{
					uint32_t k = findread_array(dict[l][ull],read,revref,revmask[j]);
					if(k != numreads)
					{
						current = k;
						flag = 1;
						updaterefcount(read[current],ref,revref,count,false,true,j);
						revcomp[i] = 1;
						sortedorder[i] = current;
						flagvec[i] = 1;//for matched
						readpos[i] = j;
						break;
					}
				}
			}	
			if(flag==1)
				break;
/*
			//look in reverse complement if nothing found yet
			std::vector<uint32_t> s1 = dict_union(revref,dict,dict_start,mask1);
			if(!s1.empty())
			{
				uint32_t k = findread(s1,read,revref,revmask[j]);
				if(k != numreads)
				{
					current = k;
					flag = 1;
					updaterefcount(read[current],ref,revref,count,false,true,j);
					revcomp[i] = 1;
					sortedorder[i] = current;
					flagvec[i] = 1;
					readpos[i] = j;
					break;
				}
			}
*/
			revref<<=2;
		}	
		
		if(flag == 0)//no match found
		{
			unmatched += 1;
			for(uint32_t j = remainingpos; ; j--)
				if(remainingreads[j] == 1)
				{
					current = j;
					remainingpos = j-1;
					break;
				}//searching from behind to speed up erase in bin
			revcomp[i] = 0;
			updaterefcount(read[current],ref,revref,count,true,false,0);
			sortedorder[i] = current;
			flagvec[i] = 0;
			readpos[i] = readlen;
		}
	}
	//remove last read and deallocate remainingreads pointer
	remainingreads[current] = 0;
	for(int l = 0; l < numdict; l++)
	{	b = read[current]&mask1[l];
		ull = (b>>2*dict_start[l]).to_ullong();
		uint32_t pos = std::lower_bound(dict[l][ull]+1,dict[l][ull]+dict[l][ull][0],current)-dict[l][ull];
		//binary search since dict[b] is sorted
		std::move(dict[l][ull]+pos+1,dict[l][ull]+dict[l][ull][0],dict[l][ull]+pos);
		dict[l][ull][0]--;
		if(dict[l][ull][0] == 1)
		{
			delete[] dict[l][ull];
			dict[l].erase(ull);
		}
	}
	delete[] remainingreads;
	std::cout << "Reordering done, "<<unmatched<<" were unmatched\n";
	return;
}


void generatemasks(std::bitset<2*readlen> *mask,std::bitset<2*readlen> *revmask)
{
	for(int i = 0; i < maxmatch; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 2*readlen - 2*i; j++)
			mask[i][j] = 1;
		for(int j = 2*i; j < 2*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}



void writetofile(std::bitset<2*readlen> *read, std::vector<uint32_t> &sortedorder,std::vector<bool> &revcomp,std::vector<bool> &flagvec, std::vector<uint32_t>& readpos)
{
	std::ofstream fout(outfile,std::ofstream::out);
	std::ofstream foutRC(outfileRC,std::ofstream::out);
	std::ofstream foutflag(outfileflag,std::ofstream::out);
	std::ofstream foutpos(outfilepos,std::ofstream::out);
	std::vector<uint32_t>::iterator it1,it4;
	std::vector<bool>::iterator it2,it3;
	char s[readlen+1],s1[readlen+1];
	s[readlen] = '\0';
	s1[readlen] = '\0';
	for (it1 = sortedorder.begin(),it2 = revcomp.begin(),it3 = flagvec.begin(),it4 = readpos.begin(); it1 != sortedorder.end(); ++it1,++it2,++it3,++it4)
	{
		foutflag << *it3;
		foutpos << *it4 << "\n";
		bitsettostring(read[*it1],s);
		if(*it2 == 0)
		{
			fout<<s<<"\n";
			foutRC << 'd';
		}
		else
		{
			reverse_complement(s,s1);
			fout<<s1<<"\n";
			foutRC << 'r';
		}
	}
	fout.close();
	foutRC.close();
	foutflag.close();
	return;
}

/*
void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	std::string s1 = b.to_string();
	char *s2 = &(s1[0]);
	for(int i = 0; i < readlen; i++)
		s[i] = inttochar[2*charinttoint[s2[2*(readlen-i-1)+1]]+charinttoint[s2[2*(readlen-i-1)]]];
	return;
}


void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	std::bitset<2*readlen> b1;
	for(int i = 0; i < readlen; i++)
	{	
		b1 = b&positionmask[i];
		if(b1 == basemask[i]['A'])
			s[i] = 'A';
		else if(b1 == basemask[i]['C'])
			s[i] = 'C';
		else if(b1 == basemask[i]['G'])
			s[i] = 'G';
		else
			s[i] = 'T';
	}
	return;
}
*/
char revinttochar[4] = {'A','G','C','T'};

void bitsettostring(std::bitset<2*readlen> b,char *s)
{
	unsigned long long ull,rem;
	for(int i = 0; i < 2*readlen/64+1; i++)
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

/*
std::bitset<2*readlen> chartobitset(char *s)
{
	char s1[2*readlen+1];
	for(int i = 0; i < readlen; i++)
	{
		s1[2*(readlen-i-1)+1] = chartobit[s[i]][0];
		s1[2*(readlen-i-1)] = chartobit[s[i]][1];
	}
	s1[2*readlen] = '\0';
	std::bitset<2*readlen> b(s1);
	return b;
}
*/

void reverse_complement(char* s, char* s1)
{
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}

void updaterefcount(std::bitset<2*readlen> cur, std::bitset<2*readlen> &ref, std::bitset<2*readlen> &revref, int count[][readlen], bool resetcount, bool rev, int shift)
{
	char s[readlen+1],s1[readlen+1],*current;
	bitsettostring(cur,s);
	if(rev == false)
		current = s;
	else
	{
		reverse_complement(s,s1);
		current = s1;
	}

	if(resetcount == true)//resetcount - unmatched read so start over
	{
		std::memset(count, 0, sizeof(count[0][0]) * 4 * readlen);	
		for(int i = 0; i < readlen; i++)
		{	
			count[chartoint[current[i]]][i] = 1;
		}

	}
	else
	{
		//shift count and find current by max
		for(int i = 0; i < readlen-shift; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = count[j][i+shift];
			count[chartoint[current[i]]][i] += 1;
			
			int max = 0,indmax = 0;
			for(int j = 0; j < 4; j++)
				if(count[j][i]>max)
				{
					max = count[j][i];
					indmax = j;
				}
			current[i] = inttochar[indmax];
		}
		//for the new positions make current same as the current read and count to 1
		for(int i = readlen-shift; i < readlen; i++)
		{	
			for(int j = 0; j < 4; j++)
				count[j][i] = 0;
			count[chartoint[current[i]]][i] = 1;
		}
	}
	ref = chartobitset(current);
	char revcurrent[readlen];
	reverse_complement(current,revcurrent);
	revref = chartobitset(revcurrent);
	return;
}
