
#ifndef KFF_IO
#define KFF_IO

// the configured options and settings for Tutorial
#define KFF_VERSION_MAJOR @KFF_VERSION_MAJOR@
#define KFF_VERSION_MINOR @KFF_VERSION_MINOR@


class Kff_file {
public:
	// encoding:       A  C  T  G
	char encoding[4] = {0, 1, 3, 2};

	Kff_file(const char * filename, const char * mode);

	// --- header functions ---
	// void set_encoding(char a, char c, char g, char t);
	// void set_metadata(size_t size, char * data);
	/* Jump to the beginning of the file and read the header.
	 * The correct encoding is set in the corresponding attribute.
	 * @return The metadata content.
	 */
	// string read_header();

	// --- general section ---
	// char read_section_name();

	// --- global variables ---
	/* Write the variable in a global variable section.
	 * If the current section is not a global var section, it close the previous section and open a gvs.
	 */
	// void write_variable(string name, uint64_t content);
	/* Read a variable from the file.
	 * The value is directly returned and the name is set in the name variable
	 */
	// uint64_t read_variable(string &s name);

	// raw section
	// minimizer section
};

#endif
