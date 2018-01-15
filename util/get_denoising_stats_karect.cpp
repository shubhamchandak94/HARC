#include <fstream>
#include <iostream>
#include <string>

std::string prefix, experiment_id;

int main(int argc, char** argv)
{
	uint64_t original_errors = 0, undetected_errors = 0, new_errors = 0, incorrectly_denoised = 0;
	prefix = std::string(argv[1]);
	experiment_id = std::string(argv[2]);

	std::ifstream f_clean(prefix+ ".upper.clean");
	std::ifstream f_noisy(prefix+ ".noisy");
	std::ifstream f_denoised(prefix+ ".karect.dna.d");
	std::string clean,noisy,denoised;
	while(std::getline(f_denoised,denoised))
	{
		std::getline(f_clean,clean);
		std::getline(f_noisy,noisy);
		for(int i = 0; i < denoised.length(); i++)
		{
			if(clean[i] != noisy[i])
			{
				original_errors++;
				if(denoised[i] == noisy[i])
					undetected_errors++;
				if((denoised[i] != noisy[i]) &&	(denoised[i] != clean[i]))
					incorrectly_denoised++;
			}
			else
				if(denoised[i] != clean[i])
					new_errors++;
		}
	}
	std::cout << "\n";
	std::cout << experiment_id << "\n";
	std::cout << "Original Errors: " << original_errors << "\n";
	std::cout << "Undetected Errors: " << undetected_errors << "\n";
	std::cout << "Incorrectly Denoised: " << incorrectly_denoised << "\n";
	std::cout << "New Errors: " << new_errors << "\n";
	std::cout << "Total Errors after Denoising; " <<  undetected_errors + incorrectly_denoised + new_errors << "\n";
	return 0;
}

