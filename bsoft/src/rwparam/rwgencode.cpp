/**
@file	rwgencode.cpp
@brief	Library routines to read and write genetic codes
@author Bernard Heymann
@date	Created: 20030316
@date	Modified: 20250601
**/

#include "rwgencode.h"
#include "star.h"
#include "mol_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
map<string,char>	read_gencode_star(string& propfile);
int 		write_genetic_code_star(string& filename, map<string,char>& gc);


/**
@brief 	Gets the genetic code from a parameter file.
@param 	&filename			file name (if empty, use a default file).
@return map<string,char>	64 element map.
**/
map<string,char>	get_genetic_code(string& filename)
{
	if ( verbose & VERB_DEBUG )	
		cout << "DEBUG get_genetic_code:" << endl;
	
	map<string,char>	code;
		
	// Atom parameter file
	string				gcfile = "gencode.star";
	if ( filename.c_str() ) gcfile = filename;
	
	string				propfile = parameter_file_path(gcfile);
	string				ext = extension(gcfile);
	if ( ext.length() ) {
		if ( ext.find("star") != string::npos )
			code = read_gencode_star(propfile);
	}
	
	if ( code.size() != 64 ) {
		if ( verbose )
			cout << "Genetic code file " << gcfile << " not opened! Using default code" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG get_genetic_code: " << code.size() << endl;
		
	return code;
}

/**
@brief 	Writing genetic code.
@param 	&filename	file name.
@param 	&gc			64 element map with the genetic code.
@return int 			0.
**/
int 		write_genetic_code(string& filename, map<string,char>& gc)
{

	return write_genetic_code_star(filename, gc);
}

int 		write_genetic_code_star(string& filename, map<string,char>& gc)
{
	string			id("Genetic_code");
	Bstar			star;

	star.comment("# Genetic code\n\n");

	BstarBlock&		block = star.add_block(id);

	BstarLoop&		loop = block.add_loop();
	loop.tags()[RESPROP_CODON] = 0;
	loop.tags()[RESPROP_CODE1] = 1;

	for ( auto c: gc ) {
		vector<string>&	vs = loop.add_row(2);
		vs[0] = c.first;
		vs[1] = c.second;
	}
	
	return star.write(filename);
}

map<string,char>	read_gencode_star(string& propfile)
{
 	Bstar				star;
	map<string,char>	gc;
	
 	if ( star.read(propfile) < 0 )
		error_show(propfile.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return gc;
	}

	int					i, j;

	for ( auto ib: star.blocks() ) {
		for ( auto il: ib.loops() ) {
			if ( ( i = il.find(RESPROP_CODE1) ) >= 0 &&
					( j = il.find(RESPROP_CODON) ) >= 0 ) {
				if ( il.data().size() != 64 ) {
					cerr <<  "Error: File " << propfile << " does contain only " << il.data().size() << " codons, 64 required!" << endl;
					return gc;
				}
				for ( auto ir: il.data() )
					gc[ir[j]] = ir[i][0];
			}
		}
	}
		
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_gencode_star: " << gc.size() << endl;
		
	return gc;
}
