/**
@file	bseq.cpp
@brief	A program to manipulate DNA and protein sequences
@author Bernard Heymann
@date	Created: 20000808
@date	Modified: 20250510
**/

#include "rwsequence.h"
#include "rwgencode.h"
#include "seq_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bseq [options] input.seq output.seq",
"------------------------------------------",
"Manipulates DNA and protein sequences.",
" ",
"Actions:",
"-show                    Show all DNA and protein sequences after other operations.",
"-Mass                    Show the calculated mass for each sequence.",
"-complement              Complement all DNA sequences before other operations.",
"-translate 1             Translate all DNA sequences with this frame (0,1,2).",
"-finddna ggattcga        DNA sequence to look for in DNA sequences.",
"-findprotein fghwereaas  Protein sequence to look for in protein sequences.",
"-findcoding aklsdrtv     Coding sequence to look for in DNA sequences.",
"                         Can be a file name with a protein sequence.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-length 1100,1500        Sequence length range to select for.",
"-number 10,14            Number of residues on either side to include in output.",
"-threshold 50            Threshold percentage to report hits (default only best hit).",
" ",
"Input:",
"-geneticcode file.star   Genetic code file (default = internal values).",
"-elements prop.star      Calculate elemental composition.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	int 			setshow(0);				// Show flag
	int 			setmass(0);				// Show flag
	int 			setfind(0);				// Search flag
	int 			setcomplement(0);		// Complement flag
	int 			settranslate(0);		// Translation flag
	int 			frame(0);				// Frame for translation
	int 			side1(0), side2(0);  	// Additional residues
	double			threshold(0);			// Only best match
	int 			seqlenmin(0);			// Minimum sequence length
	int 			seqlenmax(100000);		// Minimum sequence length
	
	// Genetic code parameter file and template sequence file
	string			gcfile;
	string			seqfile;				// Sequence to search for in a file
	string			rpfile;					// Residue properties file
	string			paramfile;
    
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "show" )
			setshow = 1;
		if ( curropt->tag == "Mass" )
			setmass = 1;
		if ( curropt->tag == "complement" )
			setcomplement = 1;
		if ( curropt->tag == "translate" ) {
        	 if ( ( frame = curropt->value.integer() ) < 0 )
				cerr << "-translate: A frame must be specified!" << endl;
			while ( frame < 0 ) frame += 3;
			while ( frame > 2 ) frame -= 3;
			settranslate = 1;
		}
		if ( curropt->tag == "length" ) {
        	if ( curropt->values(seqlenmin, seqlenmax) < 1 )
				cerr << "-length: Minimum and maximum lengths must be specified!" << endl;
			else {
				if ( seqlenmin < 0 ) seqlenmin = 0;
				if ( seqlenmax < seqlenmin ) {
					i = seqlenmin;
					seqlenmin = seqlenmax;
					seqlenmax = i;
				}
				if ( seqlenmax == seqlenmin ) seqlenmax = seqlenmin + 1;
			}
		}
		if ( curropt->tag == "finddna" ) {
			seqfile = curropt->filename().str();
			setfind = 1;
		}
		if ( curropt->tag == "findprotein" ) {
			seqfile = curropt->filename().str();
			setfind = 2;
		}
		if ( curropt->tag == "findcoding" ) {
			seqfile = curropt->filename().str();
			setfind = 3;
		}
		if ( curropt->tag == "Sequence" ) {
        	seqfile = curropt->filename().str();
			if ( seqfile.length() ) setfind = 3;
		}
		if ( curropt->tag == "number" )
			if ( curropt->values(side1, side2) < 1 )
				cerr << "-number: At least one number must be specified!" << endl;
		if ( curropt->tag == "threshold" )
        	if ( ( threshold = curropt->value.real() ) < 1e-30 )
				cerr << "-threshold: A threshold must be specified!" << endl;
		if ( curropt->tag == "geneticcode" )
            gcfile = curropt->filename().str();
		if ( curropt->tag == "elements" )
            rpfile = curropt->filename().str();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
	// Read the sequence request file if given
	vector<Bsequence>	seqs;
	string				seq;		// Sequence to search for
	if ( seqfile.length() ) {
		if ( seqfile.find(".") != string::npos ) {
			// Read the sequence file
			seqs = read_sequence(seqfile);
			if ( seqs.size() < 1 ) {
				cerr << "Error: No sequences read!" << endl;
				bexit(-1);
			}
			seq = seqs[0].sequence();
		} else {
			seq = seqfile;
		}
		// Convert threshold to number of residues
		threshold = floor(threshold*seq.length()/100.0);
	}
	
	if ( optind >= argc ) {
		cerr << "Error: No input file given!" << endl;
		bexit(-1);
	}
	
	if ( verbose && seq.length() ) cout << "Search string: " << seq << endl;
	
    // Read the sequence file
	seqfile = argv[optind++];
	seqs = read_sequence(seqfile);
	if ( seqs.size() < 1 ) {
		cerr << "Error: No sequences read!" << endl;
		bexit(-1);
	}
	
	if ( setcomplement )
		sequence_complement_all(seqs);
	
	if ( settranslate )
		sequence_translate_all(seqs, frame, gcfile);
/*		
	switch ( setfind ) {
		case 1: sequence_find_dna(seqs, seq); break;
		case 2: sequence_find_protein(seqs, seq); break;
		case 3: sequence_find_protein_in_dna(seqs, seq, seqlenmin, seqlenmax, side1, side2, threshold, gcfile); break;
		default: break;
	}
*/
	if ( setshow )
		sequence_show(seqs);
	
	if ( setmass )
		sequence_mass(seqs);
	
	if ( rpfile.length() )
		sequence_elements(seqs, rpfile);
		
	if ( optind < argc )
		write_sequence(argv[optind], seqs);
	
	timer_report(ti);
	
	bexit(0);
}

