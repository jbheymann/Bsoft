/**
@file	bcormut.cpp
@brief	A program to analyze protein sequences for correlated mutations.
@author Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250510
**/

#include "seq_analysis.h"
#include "seq_util.h"
#include "rwsequence.h"
#include "rwimg.h"
#include "rwresprop.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bcormut [options] input.pir output.pir",
"---------------------------------------------",
"Analyzes for correlated mutations in aligned protein sequences.",
"The -cutoff option must be used to do the correlated mutation analysis.",
"All sequence numbers and position numbers start at 1.",
"Sequence numbering is according to the input, including gaps.",
" ",
"Actions:",
"-cutoff 0.3              Correlated mutations: cutoff to report.",
"-limit                   Limit output to reference sequence.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype b              Set the output image data type.",
"-reference HUMAN         Reference sequence identifier (substring).",
" ",
"Input:",
"-properties file.star    Property file (default = internal values).",
" ",
"Output:",
"-matrix file.mat         Correlation coefficient matrix.",
"-image file.map          Correlation coefficient matrix as an image.",
" ",
"Examples:",
"bcormut -verbose 7 -cutoff 0.4 -image coeff.map -datatype byte seq.pir",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	string 			refseq;  				// Reference sequence ID
	int				limit(0);				// Flag to limit to reference
	double			cutoff(0.8);			// Cutoff for scoring functions
   
    // Initialize the file names and image structure
    string			matfile;
    string			imgfile;
    string			propfile;
	
	int				i, j, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "reference" ) {
        	refseq = curropt->value.str();
        	if ( refseq.length() < 1 )
				cerr << "-reference: A reference sequence ID must be specified!" << endl;
		}
		if ( curropt->tag == "limit" ) limit = 1;
		if ( curropt->tag == "cutoff" )
        	if ( ( cutoff = curropt->value.real() ) < 1e-10 )
				cerr << "-cutoff: A cutoff must be specified!" << endl;
		if ( curropt->tag == "properties" )
            propfile = curropt->filename().str();
		if ( curropt->tag == "matrix" )
			matfile = curropt->filename().str();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename().str();
    }
	option_kill(option);
    
	double		ti = timer_start();
	
    // Read the sequence file
	string		filename(argv[optind++]);
    vector<Bsequence>	seqs = read_sequence(filename);
	if ( seqs.size() < 1 ) {
		cerr << "Error: No sequences read!" << endl;
		bexit(-1);
	}
		
	long			maxlen = sequence_maximum_length(seqs);
	vector<int>		seqflag(maxlen,1);
	if ( limit ) seqflag = sequence_limit(seqs, refseq);
	
	for ( i=j=0; i<maxlen; i++ )
		j += seqflag[i];
	
	if ( verbose & VERB_PROCESS )
		cout << "Alignment positions flagged:    " << j << endl;

	Bresidue_matrix		simat;
	Matrix				mat;

	simat = get_residue_matrix(propfile);
	
	mat = sequence_correlated_mutation(seqs, seqflag, refseq, cutoff, simat);
	
	if ( optind < argc )
		write_sequence(argv[optind], seqs);
	
	if ( matfile.length() && mat.rows() )
		mat.write(matfile);

	if ( imgfile.length() && mat.rows() ) {
		Bimage*		pimg = new Bimage(mat, 1);
		pimg->change_type(nudatatype);
		write_img(imgfile, pimg, 0);
		delete pimg;
	}
	
	timer_report(ti);
	
	bexit(0);
}

