#include <iostream>
#include <cassert>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include "kff_io.hpp"
#include <math.h>
using namespace std;

void encode_sequence(char * sequence, size_t size, uint8_t * encoded);
string decode_sequence(uint8_t * encoded, size_t size);

int main(int argc, char * argv[]) {
    //arguments: -k [k-mer length] -i [input kmc file] -o [output file]
    uint64_t k = atoi(argv[1]);
    string input_path = argv[2];
    const char* output_path = argv[3];
    // --- header writing ---
	Kff_file file = Kff_file(output_path, "w");
	// Set encoding   A  C  G  T
	file.write_encoding(0, 1, 3, 2);
	// Set metadata
	file.write_metadata(11, "D@rK W@99ic");
    // --- global variable write ---

	// Set global variables
	Section_GV sgv = file.open_section_GV();
	sgv.write_var("k", k);
	sgv.write_var("max", 240);
	sgv.write_var("data_size", 1);
	sgv.close();

	// --- Write a raw sequence bloc ---
	Section_Raw sr = file.open_section_raw();
	// 2-bit sequence encoder
	uint8_t encoded[1024];
	uint8_t counts[255];
    string seq_inp;
    string count;
    std::fstream infile(input_path);
    cout << "Opened file " << input_path << "." << endl;
    while (infile >> seq_inp >> count) {
        char seq_inp_char[seq_inp.size()+1];
        strcpy(seq_inp_char, seq_inp.c_str());
        encode_sequence(seq_inp_char, k, encoded);
        counts[0] = stoi(count);
        sr.write_compacted_sequence(encoded, k, counts);
    }
    sr.close();
    /* Reading (for debugging) commented out

    file = Kff_file(output_path, "r");
	file.read_encoding();
	char metadata[uint64_t(pow(2, k))];
	uint32_t size = file.size_metadata();
	file.read_metadata(size, metadata);

	// --- Global variable read ---
	char section_name = file.read_section_type();
	sgv = file.open_section_GV();

	uint64_t max = file.global_vars["max"];
	uint64_t data_size = file.global_vars["data_size"];

	// --- Read Raw Block ---
	section_name = file.read_section_type();
	cout << "Read section " << section_name << endl;
	sr = file.open_section_raw();
	cout << "nb blocks: " << sr.nb_blocks << endl;

	uint8_t * seq = new uint8_t((max + k) / 8 + 1);
	uint8_t * data = new uint8_t(max * data_size);
	for (auto i=0 ; i<sr.nb_blocks ; i++) {
		cout << "bloc " << (i+1) << ": ";
		uint32_t nb_kmers = sr.read_compacted_sequence(seq, data);
		cout << nb_kmers << " kmers : ";
		cout << decode_sequence(seq, nb_kmers + k - 1) << " - ";
		for (auto i=0 ; i<nb_kmers ; i++)
			cout << (uint64_t)data[i] << ", ";
		cout << endl;
	}
	cout << endl;/**/
    file.close();
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
string decode_sequence(uint8_t * encoded, size_t size) {
	stringstream ss;
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

