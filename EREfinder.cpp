#include <iostream>
#include "Fasta_Loader.h"
#include <string>
#include <sstream>

int main(int argc, char* argv[])
{
	// Get arguments

	std::string query;
	std::string tempstring1, tempstring2;
	std::string infile_name, outfile_name;
	int sliding_window_width, slide_distance, max_records_parameter;
	bool save_every_ERE_value, many_output_files;
	bool interactive_mode, mask_N;
	bool alpha_equation;


	// set defaults
	infile_name = "input.fasta";
	outfile_name = "EREout.csv";
	sliding_window_width = 20000;
	slide_distance = 10000;
	max_records_parameter = 1000;
	save_every_ERE_value = false;
	many_output_files = false;
	interactive_mode = false;
	mask_N = true;
	alpha_equation = true;


	if (argc == 1)
	{
		std::cout << "\n(I)nteractive or (H)elp?\n";
		std::cin >> query;
		if (query == "H" || query == "h")
		{
			std::cout << "\nEREfinder\n";
			std::cout << "-i:\tinput file (default = input.fasta)\n";
			std::cout << "-o:\toutput file (default = EREout.csv)\n";
			std::cout << "-w:\tsliding window width (default = 20000 bp)\n";
			std::cout << "-d:\tslide distance (default = 10000 bp)\n";
			std::cout << "-r:\tmaximum number of records to run (default = 1000)\n";
			std::cout << "-v:\tsave ERE value for every bp [will be VERY large file](y or n)(default = no)\n";
			std::cout << "-m:\tmany output files [a different file for each record](y or n)(default = no)\n";
			std::cout << "-n:\tmask sequences containing N (y or n)(default = yes)\n";
			std::cout << "-a:\tuse (a)lpha or (b)eta estrogen receptor equation (default = alpha)\n";
			std::cout << "no arguments:\tinteractive mode or help\n";
			return 0;
		}
		else
		{
			// Get values for parameters for interactive mode
			std::cout << "\nEREfinder Interactive mode!";
			std::cout << "\nInput filename:\n";
			std::cin >> infile_name;

			std::cout << "Output filename:\n";
			std::cin >> outfile_name;
			
			std::cout << "Sliding window width:\n";
			std::cin >> sliding_window_width;
			
			std::cout << "Window slide interval:\n";
			std::cin >> slide_distance;

			std::cout << "Maximum number of fasta records to analyze:\n";
			std::cin >> max_records_parameter;
		
			std::cout << "Save ERE value for every base [VERY LARGE FILE!](y or n)?\n";
			std::cin >> tempstring1;
			if (!tempstring1.empty())
			{
				if(tempstring1[0] == 'y' || tempstring1[0] == 'Y')
					save_every_ERE_value = true;
			}
			std::cout << "Save a different output file for each fasta record (y or n)?\n";
			std::cin >> tempstring1;
			if (!tempstring1.empty())
			{
				if (tempstring1[0] == 'y' || tempstring1[0] == 'Y')
					many_output_files = true;
			}
			std::cout << "Mask sequences containing Ns (y or n)?\n";
			std::cin >> tempstring1;
			if (!tempstring1.empty())
			{
				if (tempstring1[0] == 'n' || tempstring1[0] == 'N')
					mask_N = false;
			}
			std::cout << "Use (a)lpha or (b)eta estrogen receptor equation?\n";
			std::cin >> tempstring1;
			if (!tempstring1.empty())
			{
				if (tempstring1[0] == 'b' || tempstring1[0] == 'B')
					alpha_equation = false;
			}


			std::cout << "\n\nParameter Values:";
			std::cout << "\nInput filename:\n" << infile_name;
			std::cout << "\nOutput filename:\n" << outfile_name;
			std::cout << "\nSliding window width:\n" << sliding_window_width;
			std::cout << "\nWindow slide interval:\n" << slide_distance;
			std::cout << "\nSave ERE value for every base (y or n)?\n";
				if (save_every_ERE_value) std::cout << "Yes";
				else std::cout << "No";
			std::cout << "\nSave a different output file for each fasta record (y or n)?\n";
				if (many_output_files) std::cout << "Yes";
				else std::cout << "No";
			std::cout << "\nMask sequences containing Ns (y or n)?\n";
				if (mask_N) std::cout << "Yes";
				else std::cout << "No";
			std::cout << "\nUse (a)lpha or (b)eta estrogen receptor equation?\n";
				if (alpha_equation) std::cout << "Alpha";
				else std::cout << "Beta";
			std::cout << "\n\nContinue (y or n)?";
			std::cin >> tempstring1;
			if (tempstring1[0] != 'y' && tempstring1[0] != 'Y')
				return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			std::cout << "\nEREfinder\n";
			std::cout << "-i:\tinput file (default = input.fasta)\n";
			std::cout << "-o:\toutput file (default = EREout.fasta)\n";
			std::cout << "-w:\tsliding window width (default = 20000 bp)\n";
			std::cout << "-d:\tslide distance (default = 10000 bp)\n";
			std::cout << "-r:\tmaximum number of records to run (default = 1000)\n";
			std::cout << "-v:\tsave ERE value for every bp [will be VERY large file](y or n)(default = no)\n";
			std::cout << "-m:\tmany output files [a different file for each record](y or n)(default = no)\n";
			std::cout << "-n:\tmask sequences containing N (y or n)(default = yes)\n";
			std::cout << "-a:\tuse (a)lpha or (b)eta estrogen receptor equation (default = alpha)\n";
			std::cout << "no arguments:\tinteractive mode or help\n";
			return 0;
		}
	}

	for (size_t i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-i")
			infile_name = tempstring2;
		if (tempstring1 == "-o")
			outfile_name = tempstring2;
		if (tempstring1 == "-w")
		{
			std::stringstream ss;
			ss << tempstring2;
			ss >> sliding_window_width;
		}
		if (tempstring1 == "-d")
		{
			std::stringstream ss;
			ss << tempstring2;
			ss >> slide_distance;
		}
		if (tempstring1 == "-r")
		{
			std::stringstream ss;
			ss << tempstring2;
			ss >> max_records_parameter;
		}
		if (tempstring1 == "-v")
		{
			if (tempstring2[0] == 'y' || tempstring2[0] == 'Y')
				save_every_ERE_value = true;
		}
		if (tempstring1 == "-m")
		{
			if (tempstring2[0] == 'y' || tempstring2[0] == 'Y')
				many_output_files = true;
		}
		if (tempstring1 == "-n")
		{
			if (tempstring2[0] == 'n' || tempstring2[0] == 'N')
				mask_N = false;
		}
		if (tempstring1 == "-a")
		{
			if (tempstring2[0] == 'b' || tempstring2[0] == 'B')
				alpha_equation = false;
		}
	}

	bool file_good;
	fastaloader my_FL;
	my_FL.max_to_load = max_records_parameter;


	file_good = my_FL.load_from_file_faster(infile_name);
	if (!file_good)
	{
		std::cout << "\nFailure to Load File!\n";
		return 0;
	}
	if (file_good)
	{
		std::cout << "\nFile Loaded!\n";
	}
	

	std::cout << "\n\nParameter Values:";
	std::cout << "\nInput filename:\n" << infile_name;
	std::cout << "\nOutput filename:\n" << outfile_name;
	std::cout << "\nSliding window width:\n" << sliding_window_width;
	std::cout << "\nWindow slide interval:\n" << slide_distance;
	std::cout << "\nSave ERE value for every base (y or n)?\n";
	if (save_every_ERE_value) std::cout << "Yes";
	else std::cout << "No";
	std::cout << "\nSave a different output file for each fasta record (y or n)?\n";
	if (many_output_files) std::cout << "Yes";
	else std::cout << "No";
	std::cout << "\nMask sequences containing Ns (y or n)?\n";
	if (mask_N) std::cout << "Yes";
	else std::cout << "No";
	std::cout << "\nUse (a)lpha or (b)eta estrogen receptor equation?\n";
	if (alpha_equation) std::cout << "Alpha\n";
	else std::cout << "Beta\n";

	// Cycle through all fasta records

	//my_FL.output_records();
	for (size_t i = 0; i < my_FL.recordlist.size(); i++)
	{
		my_FL.recordlist[i].calc_ere_vals(mask_N, alpha_equation);
		//my_FL.recordlist[i].output_ere_vals();
		my_FL.recordlist[i].calc_sliding_window(sliding_window_width, slide_distance, alpha_equation);
		//my_FL.recordlist[i].output_sliding_window();
	}

	if (save_every_ERE_value)
	{
		my_FL.save_ere_vals(outfile_name + "_ERE_values.txt");
	}

	if (!many_output_files)
		my_FL.save_sliding_window(outfile_name, alpha_equation);
	else
		my_FL.save_many_files(outfile_name, alpha_equation);

	std::cout << "\nDONE!";

	return 0;
}