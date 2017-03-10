#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <random>


void read_reference(FILE* in_file_, const long long file_size_, char* & ref_,
					long long& ref_size_, long long& no_symb_)
{
	ref_ = new char[file_size_];
	ref_size_ = 0;
	no_symb_ = 0;

	long long len = 0;

	const int row_size = 1024;
	char row[row_size];
	while (!feof(in_file_))
	{
		fgets(row, row_size, in_file_);
		len = strlen(row);

		if (len < 2)
			break;

		if (row[0] == '>')
			continue;

		memcpy(ref_ + ref_size_, row, len-1);
		ref_size_ += len-1;
	}
	ref_[ref_size_] = 0;


	// re-parse reference
	//
	for (long long i = 0; i < ref_size_; ++i)
	{
		if (ref_[i] == 'N')
			continue;

		if (ref_[i] == 'A' || ref_[i] == 'C' || ref_[i] == 'G' || ref_[i] == 'T')
			no_symb_++;
		else
			ref_[i] = 'N';
	}
}

void generate_fastq(const char* ref_, const long long ref_len_, const long long read_count_,
					const int read_len_, const bool with_errors_, FILE* out_file_)
{
	// translation lookup tables
	//
	int seq_sym[128], rc_sym[128];
	std::fill(seq_sym, seq_sym + 128, -1);
	std::fill(rc_sym, rc_sym + 128, -1);
	seq_sym['A'] = seq_sym['a'] = 0;
	seq_sym['G'] = seq_sym['g'] = 1;
	seq_sym['C'] = seq_sym['c'] = 2;
	seq_sym['T'] = seq_sym['t'] = 3;
	rc_sym['A'] = 'T';
	rc_sym['T'] = 'A';
	rc_sym['C'] = 'G';
	rc_sym['G'] = 'C';

	char seq_rand_trans[4][4] = {	{'G', 'C', 'T', 'N'},
									{'A', 'C', 'T', 'N'},
									{'A', 'G', 'T', 'N'},
									{'A', 'G', 'C', 'N'},
	};

	// FASTQ line buffers
	//
	char* read_id = new char[32];
	char* read_plus = new char[3];
	char* read_quality = new char[read_len_ + 2];
	char* read_dna = new char[read_len_ + 2];

	for (int i = 0; i < read_len_; ++i)
		read_quality[i] = 'H';
	read_quality[read_len_] = '\n';
	read_quality[read_len_+1] = 0;
	read_dna[read_len_] = '\n';
	read_dna[read_len_+1] = 0;
	read_plus[0] = '+';
	read_plus[1] = '\n';
	read_plus[2] = 0;

	std::mt19937 rnd; // Mersenne Twister random number generator


	// generate reads
	//
	for (long long i = 0; i < read_count_; )
	{
		long long ref_pos = rnd();
		if (ref_pos >= ref_len_ - read_len_)
			continue;

		if (memchr(ref_ + ref_pos, 'N', read_len_) != NULL)
			continue;

		sprintf(read_id, "@T.%lld\n", i);

		const char* ref_dna = ref_ + ref_pos;

		if (i & 1)		// try reverse-compliment ?
		{
			for	(int i = 0; i < read_len_; ++i)
				read_dna[i] = rc_sym[ref_dna[read_len_-i]];
		}
		else
		{
			for (int i = 0; i < read_len_; ++i)
				read_dna[i] = ref_dna[i];
		}

		if (with_errors_)
		{
			for (int j = 0; j < read_len_; j++)
			{
				if (rnd() % 100 != 0)
					continue;

				int c = (int)read_dna[j];
				int s = (int)(rnd() % 4);
				read_dna[j] = seq_rand_trans[seq_sym[c]][s];
			}
		}

		fputs(read_id, out_file_);
		fputs(read_dna, out_file_);
		fputs(read_plus, out_file_);
		fputs(read_quality, out_file_);

		++i;
		if (i % (4000000) == 0)
			std::cout << i / 1000000 << "M\n";
	}

	delete read_id;
	delete read_plus;
	delete read_quality;
	delete read_dna;
}


int main(int argc, char **argv)
{
	if (argc < 5)
	{
		std::cerr << "Usage: gen_fastq <no_reads> <length> <input_file> <output_file> [-e]\n";
		return 1;
	}


	// parse params
	//
	long long no_reads = atoi(argv[1]);
	int read_len = atoi(argv[2]);
	bool gen_err = (argc == 6) && strlen(argv[5]) == 2 && argv[5][1] == 'e';
	if (read_len == 0 || no_reads == 0)
	{
		std::cerr << "Error: invalid parameters";
		return 1;
	}


	// prepare IO
	//
	FILE* in_file = NULL, *out_file = NULL;
	in_file = fopen(argv[3], "rb");
	out_file = fopen(argv[4], "wb");

	if (!in_file || !out_file)
	{
		if (in_file)
			fclose(in_file);
		std::cerr << "Error: cannot open files\n";
		return 1;
	}

	setvbuf(in_file, NULL, _IOFBF, 64 << 20);
	setvbuf(out_file, NULL, _IOFBF, 64 << 20);

	fseek(in_file, 0, SEEK_END);
	long long file_size = ftell(in_file);
	fseek(in_file, 0, SEEK_SET);


	std::cout << "File size: " << file_size << "\n";
	std::cout << "Reading reference\n";

	// read reference
	//
	char* ref = NULL;
	long long ref_size = 0;
	long long no_symb = 0;
	read_reference(in_file, file_size, ref, ref_size, no_symb);

	std::cout << "No. of non-Ns: " << no_symb << "\n";


	// generate reads
	//
	std::cout << "Producing " << no_reads << " reads of length " << read_len;
	if (gen_err)
		std::cout << " with errors";
	std::cout << "\n";

	generate_fastq(ref, ref_size, no_reads, read_len, gen_err, out_file);


	// cleanup
	//
	delete ref;

	fclose(in_file);
	fclose(out_file);

	return 0;
}
