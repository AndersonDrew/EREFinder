#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

inline double convert_sequence_to_value(std::string seq) // This function not used
{
	if (seq.size() < 15)
		return -1;

	// determine perfect hits
	double N_perfect_HS = 0;
	bool perfect1 = false;
	bool perfect2 = false;
	if (seq.substr(0, 6) == "AGGTCA")
	{
		N_perfect_HS++;
		perfect1 = true;
	}
	if (seq.substr(9, 6) == "TGACCT")
	{
		N_perfect_HS++;
		perfect2 = true;
	}

	// check individual base pairs for substitutions
	double N_substitutions = 0;
	if (!perfect1)
	{
		if (seq[0] != 'A' && seq[0] != 'a' && seq[0] != 'T' && seq[0] != 't')
			N_substitutions++;
		if (seq[1] != 'G' && seq[1] != 'g' && seq[1] != 'C' && seq[1] != 'c')
			N_substitutions++;
		if (seq[2] != 'G' && seq[2] != 'g' && seq[2] != 'C' && seq[2] != 'c')
			N_substitutions++;
		if (seq[3] != 'A' && seq[3] != 'a' && seq[3] != 'T' && seq[3] != 't')
			N_substitutions++;
		if (seq[4] != 'G' && seq[4] != 'g' && seq[4] != 'C' && seq[4] != 'c')
			N_substitutions++;
		if (seq[5] != 'A' && seq[5] != 'a' && seq[5] != 'T' && seq[5] != 't')
			N_substitutions++;
	}
	if (!perfect2)
	{
		if (seq[9] != 'A' && seq[9] != 'a' && seq[9] != 'T' && seq[9] != 't')
			N_substitutions++;
		if (seq[10] != 'G' && seq[10] != 'g' && seq[10] != 'C' && seq[10] != 'c')
			N_substitutions++;
		if (seq[11] != 'A' && seq[11] != 'a' && seq[11] != 'T' && seq[11] != 't')
			N_substitutions++;
		if (seq[12] != 'G' && seq[12] != 'g' && seq[12] != 'C' && seq[12] != 'c')
			N_substitutions++;
		if (seq[13] != 'G' && seq[13] != 'g' && seq[13] != 'C' && seq[13] != 'c')
			N_substitutions++;
		if (seq[14] != 'A' && seq[14] != 'a' && seq[14] != 'T' && seq[14] != 't')
			N_substitutions++;
	}

	double ln_kd;
	ln_kd = (0.55 * N_substitutions) - (1.82 * N_perfect_HS) + 3.11;

	double kd;
	kd = exp(ln_kd);

	return kd;
} // This function not used

inline double convert_sequence_to_value_2(const std::string &seq, int start_pos, bool mask_n, bool alpha)
{
	if (start_pos + 15 > seq.size())
		return -1;

	// determine perfect hits
	double N_perfect_HS = 0;
	bool perfect1 = false;
	bool perfect2 = false;

	bool has_Ns;
	if (mask_n)
	{
		has_Ns = false;
		for (size_t i = 0; i < 6; i++)
		{
			if (seq[start_pos + i] == 'N' || seq[start_pos + i + 9] == 'N' || seq[start_pos + i] == 'n' || seq[start_pos + i + 9] == 'n')
				has_Ns = true;
		}
		if (has_Ns)
			return -1;
	}

	if ((seq[start_pos] == 'A' || seq[start_pos] == 'a') &&
		(seq[start_pos + 1] == 'G' || seq[start_pos + 1] == 'g') &&
		(seq[start_pos + 2] == 'G' || seq[start_pos + 2] == 'g') &&
		(seq[start_pos + 3] == 'T' || seq[start_pos + 3] == 't') &&
		(seq[start_pos + 4] == 'C' || seq[start_pos + 4] == 'c') &&
		(seq[start_pos + 5] == 'A' || seq[start_pos + 5] == 'a'))
	{
		N_perfect_HS++;
		perfect1 = true;
	}
	//if (seq.substr(start_pos + 9, 6) == "TGACCT")
	if ((seq[start_pos+9] == 'T' || seq[start_pos + 9] == 't') &&
		(seq[start_pos + 10] == 'G' || seq[start_pos + 10] == 'g') &&
		(seq[start_pos + 11] == 'A' || seq[start_pos + 11] == 'a') &&
		(seq[start_pos + 12] == 'C' || seq[start_pos + 12] == 'c') &&
		(seq[start_pos + 13] == 'C' || seq[start_pos + 13] == 'c') &&
		(seq[start_pos + 14] == 'T' || seq[start_pos + 14] == 't'))
	{
		N_perfect_HS++;
		perfect2 = true;
	}

	// check individual base pairs for substitutions
	double N_substitutions = 0;
	if (!perfect1)
	{
		if (seq[start_pos] != 'A' && seq[start_pos] != 'a' && seq[start_pos] != 'T' && seq[start_pos] != 't')
			N_substitutions++;
		if (seq[start_pos + 1] != 'G' && seq[start_pos + 1] != 'g' && seq[start_pos + 1] != 'C' && seq[start_pos + 1] != 'c')
			N_substitutions++;
		if (seq[start_pos + 2] != 'G' && seq[start_pos + 2] != 'g' && seq[start_pos + 2] != 'C' && seq[start_pos + 2] != 'c')
			N_substitutions++;
		if (seq[start_pos + 3] != 'A' && seq[start_pos + 3] != 'a' && seq[start_pos + 3] != 'T' && seq[start_pos + 3] != 't')
			N_substitutions++;
		if (seq[start_pos + 4] != 'G' && seq[start_pos + 4] != 'g' && seq[start_pos + 4] != 'C' && seq[start_pos + 4] != 'c')
			N_substitutions++;
		if (seq[start_pos + 5] != 'A' && seq[start_pos + 5] != 'a' && seq[start_pos + 5] != 'T' && seq[start_pos + 5] != 't')
			N_substitutions++;
	}
	if (!perfect2)
	{
		if (seq[start_pos + 9] != 'A' && seq[start_pos + 9] != 'a' && seq[start_pos + 9] != 'T' && seq[start_pos + 9] != 't')
			N_substitutions++;
		if (seq[start_pos + 10] != 'G' && seq[start_pos + 10] != 'g' && seq[start_pos + 10] != 'C' && seq[start_pos + 10] != 'c')
			N_substitutions++;
		if (seq[start_pos + 11] != 'A' && seq[start_pos + 11] != 'a' && seq[start_pos + 11] != 'T' && seq[start_pos + 11] != 't')
			N_substitutions++;
		if (seq[start_pos + 12] != 'G' && seq[start_pos + 12] != 'g' && seq[start_pos + 12] != 'C' && seq[start_pos + 12] != 'c')
			N_substitutions++;
		if (seq[start_pos + 13] != 'G' && seq[start_pos + 13] != 'g' && seq[start_pos + 13] != 'C' && seq[start_pos + 13] != 'c')
			N_substitutions++;
		if (seq[start_pos + 14] != 'A' && seq[start_pos + 14] != 'a' && seq[start_pos + 14] != 'T' && seq[start_pos + 14] != 't')
			N_substitutions++;
	}

	double ln_kd;

	if (alpha)
		ln_kd = (0.55 * N_substitutions) - (1.82 * N_perfect_HS) + 3.11;
	else
		ln_kd = (0.50 * N_substitutions) - (1.48 * N_perfect_HS) + 3.41;

	double kd;
	kd = exp(ln_kd);

	return kd;
}


class fastarecord
{
public:
	std::string header;
	std::string sequence;
	std::vector<int> sw_sample_size;
	std::vector<double> ere_values;
	std::vector<int> sw_position;
	std::vector<double> sw_value;
	std::vector<int> N_perfect;
	std::vector<int> N_1_mismatch;
	std::vector<int> N_2_mismatch;
	std::vector<int> N_3_mismatch;
	std::vector<int> N_total;

	void calc_ere_vals(bool mask_n, bool alpha)
	{
		int i;
		double temp_ere;
		ere_values.clear();
		ere_values.reserve(sequence.size());
		for (i = 0; i < sequence.size() - 14; i++)
		{
			temp_ere = convert_sequence_to_value_2(sequence, i, mask_n, alpha);
			ere_values.push_back(1/temp_ere);
			if ((i + 1) % 1000000 == 0)
				std::cout << "\nBasepairs done:\t" << i + 1;
		}
	}

	void output_ere_vals()
	{
		std::cout << "\n\n";
		for (size_t i = 0; i < ere_values.size(); i++)
		{
			std::cout << ere_values[i] << ",";
		}
	}

	void calc_sliding_window(int window_size, int slide_distance, bool alpha)
	{
		bool off_the_end;
		off_the_end = false;
		int slide_counter = 0;
		int start, end;
		sw_sample_size.clear();
		sw_position.clear();
		sw_value.clear();
		N_perfect.clear();
		N_1_mismatch.clear();
		N_2_mismatch.clear();
		N_3_mismatch.clear();
		N_total.clear();

		while (!off_the_end)
		{
			start = slide_counter*slide_distance;
			end = slide_counter*slide_distance + window_size;
			if (end > static_cast<int>(ere_values.size()))
			{
				end = static_cast<int>(ere_values.size());
				off_the_end = true;
			}

			if (start < end)
			{
				double mean = 0;
				double number_in_window = 0;
				int perf_count = 0;
				int N_1_count = 0;
				int N_2_count = 0;
				int N_3_count = 0;
				for (int i = start; i < end; i++)
				{
					if (ere_values[i] > -1)
					{
						mean = mean + ere_values[i];
						number_in_window++;
						if (alpha)
						{
							if (ere_values[i] >= 1.5) perf_count++;
							if (ere_values[i] >= 0.15 && ere_values[i] < 1.5) N_1_count++;
							if (ere_values[i] >= 0.09 && ere_values[i] < 0.15) N_2_count++;
							if (ere_values[i] >= 0.04 && ere_values[i] < 0.09) N_3_count++;
						}
						else
						{
							if (ere_values[i] >= 0.5) perf_count++;
							if (ere_values[i] >= 0.08 && ere_values[i] < 0.5) N_1_count++;
							if (ere_values[i] >= 0.05 && ere_values[i] < 0.08) N_2_count++;
							if (ere_values[i] >= 0.03 && ere_values[i] < 0.05) N_3_count++;
						}
					}
				}
				if (number_in_window > 0)
					mean = mean / number_in_window;
				else
					mean = -1;

				sw_sample_size.push_back(static_cast<int>(number_in_window));
				sw_position.push_back(start + window_size / 2);
				sw_value.push_back(mean);
				N_perfect.push_back(perf_count);
				N_1_mismatch.push_back(N_1_count);
				N_2_mismatch.push_back(N_2_count);
				N_3_mismatch.push_back(N_3_count);
				N_total.push_back(perf_count + N_1_count + N_2_count + N_3_count);
			}
			slide_counter++;
		}
	}

	void output_sliding_window()
	{
		std::cout << "\n\nSliding Window:";
		for (size_t i = 0; i < sw_position.size(); i++)
		{
			std::cout << "\n" << sw_position[i] << "\t" << sw_value[i];
			std::cout << "\t" << N_perfect[i] << "\t" << N_1_mismatch[i] << "\t";
			std::cout << N_2_mismatch[i] << "\t" << N_3_mismatch[i] << "\t" << N_total[i];
		}
		std::cout << "\n";
	}

};


class fastaloader
{
public:
	std::vector<fastarecord> recordlist;
	unsigned int max_to_load;

	bool load_from_file(std::string filename)
	{
		fastarecord temp_rec;
		bool good_load = true;
		recordlist.clear();

		//convert filename to null terminated char array
		char filename_char[256];
		size_t k;
		for (k = 0; k < filename.size(); k++)
		{
			filename_char[k] = filename[k];
		}
		filename_char[k] = '\0';

		std::ifstream infile;
		infile.open(filename_char);
		if (!infile.good())
			return false;
		
		char c;

		bool in_header;
		while (infile.get(c) && recordlist.size() < max_to_load)
		{
			if (c == '>')
			{
				in_header = true;
				if (temp_rec.header.size() > 0 && temp_rec.sequence.size() > 0) // make sure there's a record (won't run at first >)
				{
					recordlist.push_back(temp_rec);
					temp_rec.header.clear();
					temp_rec.sequence.clear();
				}
			}
			if (in_header && (c == '\r' || c == '\n'))
			{
				in_header = false;
			}

			if (in_header)
			{
				if (static_cast<int>(c) > 31)
					temp_rec.header = temp_rec.header + c;
			}
			else
			{
				if (static_cast<int>(c) > 31)
				{
					temp_rec.sequence = temp_rec.sequence + c;
					if (temp_rec.sequence.size() % 1000000 == 0)
						std::cout << "\nanother million base pairs...";
				}
			}
			//std::cout << c;
		}
		if (infile.eof()) // if it read all the way to the end, then the last record was not added to the list
		{
			recordlist.push_back(temp_rec);
		}

		//std::cout << "\nSize_reclist:\t" << recordlist.size() << "\n";

		infile.close();

		return good_load;
	}

	bool load_from_file_faster(std::string filename)
	{
		int h_counter, s_counter;

		fastarecord temp_rec;
		bool good_load = true;
		recordlist.clear();

		//convert filename to null terminated char array
		char filename_char[256];
		size_t k;
		for (k = 0; k < filename.size(); k++)
		{
			filename_char[k] = filename[k];
		}
		filename_char[k] = '\0';

		std::ifstream infile;
		infile.open(filename_char);
		if (!infile.good())
			return false;

		// Determine the size of the file and load the whole thing into a character array
		infile.seekg(0, infile.end);
		int file_length = static_cast<int>(infile.tellg());
		infile.seekg(0, infile.beg); 

		std::cout << "\nFile Size:\t" << file_length;

		char *file_buffer = new char[file_length + 1];
		char *temp_seq = new char[file_length + 1];
		char *temp_header = new char[10000];

		std::cout << "\nCharacter Array Allocated";
		
		infile.read(file_buffer, file_length);

		std::cout << "\nFile Loaded";

		int preview_count = 500;
		if (preview_count > file_length)
			preview_count = file_length;
		std::cout << "\nFirst " << preview_count << "characters of file:\n";
		for (int i = 0; i < preview_count; i++)
			std::cout << file_buffer[i];

		int buff_position;
		bool in_header;

		buff_position = 0;
		while (file_buffer[buff_position] != '>' && buff_position < file_length)
		{
			std::cout << file_buffer[buff_position];
			buff_position++;
			std::cout << "\nLooking for first > ... (check file)";
		}

		h_counter = 0;
		s_counter = 0;
		while (buff_position < file_length && recordlist.size() < max_to_load)
		{
			if (file_buffer[buff_position] == '>')
			{
				in_header = true;
				if (h_counter > 0 && s_counter > 0) // make sure there's a record (won't run at first >)
				{
					temp_header[h_counter] = '\0';
					temp_seq[s_counter] = '\0';
					temp_rec.header.assign(temp_header);
					temp_rec.sequence.assign(temp_seq);
					recordlist.push_back(temp_rec);
					temp_rec.header.clear();
					temp_rec.sequence.clear();
					h_counter = 0;
					s_counter = 0;
				}
			}
			if (in_header && (file_buffer[buff_position] == '\r' || file_buffer[buff_position] == '\n'))
			{
				in_header = false;
			}

			if (in_header)
			{
				if (static_cast<int>(file_buffer[buff_position]) > 31 && h_counter < 9999)
				{
					temp_header[h_counter] = file_buffer[buff_position];
					h_counter++;
				}
			}
			else
			{
				if (static_cast<int>(file_buffer[buff_position]) > 31)
				{
					temp_seq[s_counter] = file_buffer[buff_position];
					s_counter++;
					if (s_counter % 1000000 == 0)
						std::cout << "\nanother million base pairs...";
				}
			}
			buff_position++;
		}
		if (buff_position == file_length) // if it read all the way to the end, then the last record was not added to the list
		{
			temp_header[h_counter] = '\0';
			temp_seq[s_counter] = '\0';
			temp_rec.header.assign(temp_header);
			temp_rec.sequence.assign(temp_seq);
			recordlist.push_back(temp_rec);
			temp_rec.header.clear();
			temp_rec.sequence.clear();
			h_counter = 0;
			s_counter = 0;
		}

		//std::cout << "\nSize_reclist:\t" << recordlist.size() << "\n";

		infile.close();
		delete[] file_buffer;
		delete[] temp_seq;
		delete[] temp_header;

		return good_load;
	} // end of load_from_file_faster

	void output_records()
	{
		for (size_t i = 0; i < recordlist.size(); i++)
		{
			std::cout << recordlist[i].header << "\n";
			std::cout << recordlist[i].sequence << "\n";
		}
	}
	
	void save_ere_vals(std::string outfile)
	{
		// convert outfile to null terminated char*
		char filename[256];
		size_t i;
		for (i = 0; i < outfile.size(); i++)
			filename[i] = outfile[i];
		filename[i] = '\0';

		std::ofstream out;
		out.open(filename);
		for (i = 0; i < recordlist.size(); i++)
		{
			out << recordlist[i].header << "\n";
			for (size_t j = 0; j < recordlist[i].ere_values.size(); j++)
			{
				out << recordlist[i].ere_values[j] << ",";
			}
			out << "\n";
		}
		out.close();

	}

	void save_sliding_window(std::string outfile, bool alpha)
	{
		// convert outfile to null terminated char*
		char filename[256];
		size_t i, j;
		for (i = 0; i < outfile.size(); i++)
			filename[i] = outfile[i];
		filename[i] = '\0';

		std::ofstream out;
		out.open(filename);
		for (j = 0; j < recordlist.size(); j++)
		{
			out << recordlist[j].header << "\n";
			if (alpha)
				out << "genome_position,N,mean_Kd_inverse,N_over_1.5,N_0.15_to_1.5,N_0.09_to_0.15,N_0.04_to_0.09,Total_over_0.04";
			else
				out << "genome_position,N,mean_Kd_inverse,N_over_0.5,N_0.08_to_0.5,N_0.05_to_0.08,N_0.03_to_0.05,Total_over_0.03";
			for (i = 0; i < recordlist[j].sw_position.size(); i++)
			{
				out << "\n" << recordlist[j].sw_position[i] << "," << recordlist[j].sw_sample_size[i] << "," << recordlist[j].sw_value[i]
				<< "," << recordlist[j].N_perfect[i] << "," << recordlist[j].N_1_mismatch[i] << ","
				<< recordlist[j].N_2_mismatch[i] << "," << recordlist[j].N_3_mismatch[i] << "," << recordlist[j].N_total[i];
			}
			out << "\n";
		}
		out.close();
	}

	void save_many_files(std::string outfile, bool alpha)
	{
		size_t i, j;
		std::stringstream ss;
		std::string filename_long;
		for (j = 0; j < recordlist.size(); j++)
		{
			ss.clear();
			ss.str("");
			ss << outfile << "_fastarec_" << j + 1 << ".csv";
			ss >> filename_long;

			// convert filename_long to null terminated char*
			char filename[256];
			for (i = 0; i < filename_long.size(); i++)
				filename[i] = filename_long[i];
			filename[i] = '\0';

			std::ofstream out;
			out.open(filename);
			out << recordlist[j].header << "\n";
			if (alpha)
				out << "genome_position,N,mean_Kd_inverse,N_over_1.5,N_0.15_to_1.5,N_0.09_to_0.15,N_0.04_to_0.09,Total_over_0.04";
			else
				out << "genome_position,N,mean_Kd_inverse,N_over_0.5,N_0.08_to_0.5,N_0.05_to_0.08,N_0.03_to_0.05,Total_over_0.03";
			for (i = 0; i < recordlist[j].sw_position.size(); i++)
			{
				out << "\n" << recordlist[j].sw_position[i] << "," << recordlist[j].sw_sample_size[i] << "," << recordlist[j].sw_value[i]
					<< "," << recordlist[j].N_perfect[i] << "," << recordlist[j].N_1_mismatch[i] << ","
					<< recordlist[j].N_2_mismatch[i] << "," << recordlist[j].N_3_mismatch[i] << "," << recordlist[j].N_total[i];
			}
			out << "\n";
			out.close();
		}
		
	}
	
};