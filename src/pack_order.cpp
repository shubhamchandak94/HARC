#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>

std::string infile;

void pack_order();//pack order into least number of bits possible

int main(int argc, char** argv)
{
	std::string basedir = std::string(argv[1]);
	infile = basedir + "/read_order.bin";
	pack_order();
	return 0;
}

void pack_order()
{
	std::ofstream f_out(infile+".tmp",std::ios::binary);
	std::ifstream f_in(infile,std::ios::binary);
	std::ofstream f_tail(infile+".tail",std::ios::binary);
	
	uint32_t numreads = 0;
	uint32_t order;
	f_in.read((char*)&order, sizeof(uint32_t));
	while(!f_in.eof())
	{
		numreads++;
		f_in.read((char*)&order, sizeof(uint32_t));
	}
	f_in.close();
	f_in.open(infile,std::ios::binary);
	int numbits = (int)(log2(numreads)+1);
	f_out.write((char*)&numbits, sizeof(int));
	f_out.write((char*)&numreads, sizeof(uint32_t));
	uint32_t *order_array = new uint32_t [numbits];
	for(uint32_t i = 0; i < numreads/32; i++)
	{
		std::fill(order_array,order_array+numbits,0);
		int pos_in_array = 0, pos_in_int = 0;
		for(int k = 0; k < 32; k++)
		{
			f_in.read((char*)&order, sizeof(uint32_t));

			order_array[pos_in_array] |= (order << pos_in_int);
			if(pos_in_int+numbits > 32)
			{
				pos_in_array++;
				order_array[pos_in_array] = (order >> (32 - pos_in_int));
				pos_in_int = numbits - (32 - pos_in_int);
			}
			else if(pos_in_int+numbits == 32)
			{
				pos_in_array++;
				pos_in_int = 0;
			}
			else
				pos_in_int = pos_in_int + numbits;
		}
		for(int k = 0; k < numbits; k++)
			f_out.write((char*)&order_array[k], sizeof(uint32_t));
	}
	for(uint32_t i = 0; i < numreads%32; i++)
	{
		f_in.read((char*)&order, sizeof(uint32_t));
		f_tail.write((char*)&order, sizeof(uint32_t));	
	}
	f_tail.close();
	f_in.close();
	f_out.close();
	remove(infile.c_str());
	rename((infile+".tmp").c_str(),infile.c_str());
	return;
}
