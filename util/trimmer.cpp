#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

int main(int argc, char **argv)
{
	std::string infile_1 = std::string(argv[1]);
	std::string infile_2 = std::string(argv[2]);
	int max_readlen = atoi(argv[3]);
	std::string s[2][4];
	std::string outfile_1 = infile_1 + ".trimmed";	
	std::string outfile_2 = infile_2 + ".trimmed";
	std::ifstream f_in_1(infile_1);
	std::ifstream f_in_2(infile_2);
	std::ofstream f_out_1(outfile_1);
	std::ofstream f_out_2(outfile_2);
	while(std::getline(f_in_1,s[0][0]))
	{
		std::getline(f_in_1,s[0][1]);
		std::getline(f_in_1,s[0][2]);
		std::getline(f_in_1,s[0][3]);
		std::getline(f_in_2,s[1][0]);
		std::getline(f_in_2,s[1][1]);
		std::getline(f_in_2,s[1][2]);
		std::getline(f_in_2,s[1][3]);
		if(s[0][1].length() >= max_readlen && s[1][1].length() >= max_readlen)
		{
			f_out_1 << s[0][0] <<"\n" <<  s[0][1].substr(0,max_readlen) <<"\n" <<  s[0][2] <<"\n" <<  s[0][3].substr(0,max_readlen) <<"\n"; 
			f_out_2 << s[1][0] <<"\n" <<  s[1][1].substr(0,max_readlen) <<"\n" <<  s[1][2] <<"\n" <<  s[1][3].substr(0,max_readlen) <<"\n"; 

		}
	}
	f_in_1.close();
	f_in_2.close();
	f_out_1.close();
	f_out_2.close();
	return 0;	
}	
