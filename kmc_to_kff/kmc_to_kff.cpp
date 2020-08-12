/*
  This file is a part of KMC software distributed under GNU GPL 3 licence.
  The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

  This file demonstrates the example usage of kmc_api software. 
  It reads kmer_counter's output and prints kmers to an output file.

  Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

  Version: 3.1.1
  Date   : 2019-05-19
*/
#include <iostream>
#include "kmc_api/kmc_file.h"
#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include "kff_io.hpp"
#include <math.h>

void print_info(void);
void encode_sequence(char * sequence, size_t size, uint8_t * encoded);
std::string decode_sequence(uint8_t * encoded, size_t size);

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

int _tmain(int argc, char* argv[])
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		print_info();
		return 0;
	}

	CKMCFile kmer_data_base;
	int32 i;
	uint32 min_count_to_set = 0;
	uint32 max_count_to_set = 0;
	std::string input_file_name;
	std::string output_file_name;

	FILE * out_file;
	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 3)
	{
		print_info();
		return EXIT_FAILURE;
	}

	for(i = 1; i < argc; ++i)
	{
		if(argv[i][0] == '-')
		{	
			if(strncmp(argv[i], "-ci", 3) == 0)
				min_count_to_set = atoi(&argv[i][3]);
			else if(strncmp(argv[i], "-cx", 3) == 0)
					max_count_to_set = atoi(&argv[i][3]);
		}
		else
			break;
	}

	if(argc - i < 2)
	{ 
		print_info();
		return EXIT_FAILURE;
	}

	input_file_name = std::string(argv[i++]);
	output_file_name = std::string(argv[i]);

	if((out_file = fopen (output_file_name.c_str(),"wb")) == NULL)
	{
		print_info();
		return EXIT_FAILURE;
	}

	setvbuf(out_file, NULL ,_IOFBF, 1 << 24);

	//------------------------------------------------------------------------------
	// Open kmer database for listing and print kmers within min_count and max_count
	//------------------------------------------------------------------------------

	if (!kmer_data_base.OpenForListing(input_file_name))
	{
		print_info();
		return EXIT_FAILURE ;
	}
	else
	{
		uint32 _kmer_length;
		uint32 _mode;
		uint32 _counter_size;
		uint32 _lut_prefix_length;
		uint32 _signature_len;
		uint32 _min_count;
		uint64 _max_count;
		uint64 _total_kmers;

		kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

		Kff_file file = Kff_file(output_file_name.c_str(), "w");
	    // Set encoding   A  C  G  T
	    file.write_encoding(0, 1, 3, 2);
	    // Set metadata
	    file.write_metadata(11, "D@rK W@99ic");
        // --- global variable write ---

	    // Set global variables
	    Section_GV sgv = file.open_section_GV();
        sgv.write_var("k", _kmer_length);
        sgv.write_var("max", 240);
        sgv.write_var("data_size", 1);
        sgv.close();

	// --- Write a raw sequence bloc ---
	    Section_Raw sr = file.open_section_raw();
	// 2-bit sequence encoder
        uint8_t encoded[1024];
        uint8_t counts[255];
        sr.close();
		//std::string str;
		char str[1024];
		uint32 counter_len;
		
		CKmerAPI kmer_object(_kmer_length);

		if (_mode) //quake compatible mode
		{
			float counter;
			while (kmer_data_base.ReadNextKmer(kmer_object, counter))
			{
				kmer_object.to_string(str);
                encode_sequence(str, _kmer_length, encoded);
                counts[0] = int(counter);
                sr.write_compacted_sequence(encoded, _kmer_length, counts);		
			}
		}
		else
		{
			uint64 counter;
			while (kmer_data_base.ReadNextKmer(kmer_object, counter))
			{
				kmer_object.to_string(str);
                encode_sequence(str, _kmer_length, encoded);
                counts[0] = int(counter);
                sr.write_compacted_sequence(encoded, _kmer_length, counts);	
			}
		}
		
	
		fclose(out_file);
		kmer_data_base.Close();
	}

	return EXIT_SUCCESS; 
}
// -------------------------------------------------------------------------
// Print execution options 
// -------------------------------------------------------------------------
void print_info(void)
{
	std::cout << "KMC2KFF" << KMC_VER << " (" << KMC_DATE << ")\n"
			  << "\nUsage:\nkmc2kff[options] <kmc_database> <output_file>\n"
			  << "Parameters:\n"
			  << "<kmc_database> - kmer_counter's output\n";
}
// ---- DNA encoding functions [ACTG] -----

uint8_t uint8_packing(char * sequence, size_t size);
/* Encode the sequence into an array of uint8_t packed sequence slices.
 * The encoded sequences are organised in big endian order.
 */
void encode_sequence(char * sequence, size_t size, uint8_t * encoded) {
	// Encode the truncated first 8 bits sequence
	size_t remnant = size % 4;
	if (remnant > 0) {
		encoded[0] = uint8_packing(sequence, remnant);
		encoded += 1;
	}

	// Encode all the 8 bits packed
	size_t nb_uint_needed = size / 4;
	for (size_t i=0 ; i<nb_uint_needed ; i++) {
		encoded[i] = uint8_packing(sequence + remnant + (i<<2), 4);
	}
}

/* Transform a char * sequence into a uint8_t 2-bits/nucl
 * Encoding ACTG
 * Size must be <= 4
 */
uint8_t uint8_packing(char * sequence, size_t size) {
	assert(size <= 4);

	uint8_t val = 0;
	for (size_t i=0 ; i<size ; i++) {
		val <<= 2;
		val += (sequence[i] >> 1) & 0b11;
	}

	return val;
}


void uint8_unpacking(uint8_t packed, char * decoded, size_t size);
std::string decode_sequence(uint8_t * encoded, size_t size) {
	std::stringstream ss;
	char tmp_chars[4];

	// Decode the truncated first compacted 8 bits
	size_t remnant = size % 4;
	if (remnant > 0) {
		uint8_unpacking(encoded[0], tmp_chars, remnant);
		for (size_t i=0 ; i<remnant ; i++) {
			ss << tmp_chars[i];
		}
		encoded += 1;
	}

	// Decode all the 8 bits packed
	size_t nb_uint_used = size / 4;
	for (size_t i=0 ; i<nb_uint_used ; i++) {
		uint8_unpacking(encoded[i], tmp_chars, 4);
		for (size_t i=0 ; i<4 ; i++) {
			ss << tmp_chars[i];
		}
	}

	return ss.str();
}

char const_nucleotides[4] = {'A', 'C', 'T', 'G'};
void uint8_unpacking(uint8_t packed, char * decoded, size_t size) {
	assert(size <= 4);

	size_t offset = 4 - size;
	for (size_t i=0 ; i<size ; i++) {
		decoded[i] = const_nucleotides[(packed >> ((3-i-offset) * 2)) & 0b11];
	}
}

// ***** EOF
