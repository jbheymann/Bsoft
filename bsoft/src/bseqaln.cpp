/**
@file	bseqaln.cpp
@brief	A program to analyze aligned protein sequences
@author Bernard Heymann
@date	Created: 19990123
@date 	Modified: 20250510
**/

#include "rwsequence.h"
#include "rwresprop.h"
#include "seq_align.h"
#include "seq_analysis.h"
#include "seq_util.h"
#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bseqaln [options] input.seq output.seq",
"---------------------------------------------",
"Analyzes aligned protein sequences.",
"All sequence numbers and position numbers start at 1.",
"Numbering is according to the input, including gaps.",
" ",
"Actions:",
"-pair                    Align a pair of sequences (first two in input file).",
"-nogaps                  Strip gaps from input sequences.",
"-consolidategaps         Remove positions in alignment with only gaps.",
"-identity                Pairwise identity.",
"-similarity 0.5          Pairwise similarity: threshold for counting similar residues (default 0.5).",
"-profile                 Generate a profile from an alignment.",
"-information             Information analysis.",
"-hydrophobicity 0.8      Hydrophobicity analysis: representation threshold (fraction of # of seqs).",
" ",
"Selections:",
"-select 3,0.85           Select relative to a reference sequence and cutoff.",
"-length 128,435          Select sequences within a length range.",
"-delete                  Delete non-selected sequences.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-gap 10,0.2              Gap opening and extending penalties for alignment (default 20,0.2).",
"-segments 1-33,36-69     Selection of sequence segments (default all).",
"                         Every pair of numbers indicate a segment to select.",
"-sequences 2,18          Subset of sequences to compare (default all).",
"-window 17               Window for moving averages (default 20).",
" ",
"Input:",
"-properties file.star    Property file (default = internal values).",
"-Similarity sim.star     Residue similarity matrix (default blosum62.star).",
" ",
"Output:",
"-Matrix file.mat         Matrix output filename.",
"-Image file.img          Image output filename.",
"-Postscript file.ps      Postscript output filename.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    // Initialize variables
	int 			sub1(0), sub2(0);  		// Subset to compare with
	bool			consolidate_gaps(0);
	bool			do_pair(0);
	bool			do_id(0);
	bool			do_sim(0);
	bool			do_profile(0);
	bool			do_hp(0);
	bool			do_info(0);
	bool			delseq(0);				// Flag to delete non-selected sequences
	long			ref(0);					// Reference sequence for selection
	double			cutoff(-100);			// Selection threshold
	long			len1(0), len2(0);		// Length range to select for
	double			gapopen(20), gapextend(0.2);	// Gap opening and extending penalties
//	int 			sel[20];				// Residue range selection
	int 			window(20);				// Window for moving average
	double			sim_threshold(0.5);		// Similarity threshold
	double			rep_threshold(0);		// Representation threshold
	string			segments;				// Position selector
    
    // Initialize the file names and image structure
    string			matfile;
    string			imgfile;
    string			psfile;
	string			propfile;
	string			simfile;				// Use default similarity matrix file
	
	int				i, j, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "pair" )
        	do_pair = 1;
		if ( curropt->tag == "consolidategaps" )
        	consolidate_gaps = 1;
		if ( curropt->tag == "identity" )
        	do_id = 1;
		if ( curropt->tag == "similarity" ) {
			if ( ( sim_threshold = curropt->value.real() ) < 0.01 )
				cerr << "-similarity: A threshold must be specified!" << endl;
        	else
				do_sim = 1;
		}
		if ( curropt->tag == "profile" )
        	do_profile = 1;
		if ( curropt->tag == "information" )
        	do_info = 1;
		if ( curropt->tag == "hydrophobicity" ) {
			if ( ( rep_threshold = curropt->value.real() ) < 0.001 )
				cerr << "-hydrophobicity: A threshold must be specified!" << endl;
        	else
				do_hp = 1;
		}
		if ( curropt->tag == "window" )
			if ( ( window = curropt->value.integer() ) < 1 )
				cerr << "-window: A window size must be specified!" << endl;
		if ( curropt->tag == "gap" )
			if ( curropt->values(gapopen, gapextend) < 1 )
				cerr << "-gap: At least one gap penalty must be specified!" << endl;
		if ( curropt->tag == "segments" ) {
			segments = curropt->value.str();
			if ( segments.length() < 1 )
				cerr << "-segments: Numbers must be specified!" << endl;
		}
		if ( curropt->tag == "sequences" ) {
			if ( curropt->values(sub1, sub2) < 2 )
				cerr << "-sequences: Two numbers must be specified!" << endl;
			else {
				if ( sub1 > 0 ) sub1--;
				if ( sub2 > sub1 ) sub2--;
				else sub2 = sub1;
			}
		}
		if ( curropt->tag == "delete" ) delseq = 1;
		if ( curropt->tag == "select" )
			if ( curropt->values(ref, cutoff) < 2 )
				cerr << "-select: Two numbers must be specified!" << endl;
		if ( curropt->tag == "length" )
			if ( curropt->values(len1, len2) < 2 )
				cerr << "-length: Two lengths must be specified!" << endl;
		if ( curropt->tag == "properties" )
			propfile = curropt->filename().str();
		if ( curropt->tag == "Similarity" )
			simfile = curropt->filename().str();
		if ( curropt->tag == "Matrix" )
			matfile = curropt->filename().str();
		if ( curropt->tag == "Image" )
			imgfile = curropt->filename().str();
		if ( curropt->tag == "Postscript" )
			psfile = curropt->filename().str();
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
	
//	cout << "Properties file = " << propfile << endl;
	if ( do_hp ) propfile.clear();
	
//	molecule_update_comment(molgroup, argc, argv);
	
//	Bmolecule*	mol = NULL;
	string		profile;

	// Check the moving average window size
	if ( window < 2 ) window = 2;
	if ( window > 1000 ) window = 1000;
	
	// Set the tags according to the sequences selected
	long			maxlen = sequence_maximum_length(seqs);
	vector<int>		seqflag(maxlen,1);
	if ( sub2 > 0 ) {
		for ( j=0; j<maxlen; j++ )
			seqflag[j] = 0;
		i = 0;
		for ( auto seq: seqs ) {
			if ( i>=sub1 && i<=sub2 ) {
				string& 	seqstr = seq.sequence();
				for ( j=0; j<seqstr.length(); j++ )
					if ( seqstr[j] != '-' )
						seqflag[j] = 1;
			}
			i++;
		}
	}
	
	// Get the selection
/*	if ( segments.length() ) {
		pnt = segments;
		n = 0;
		while ( pnt && sscanf(pnt, "%d", &sel[n]) ) {
			if ( ( pnt = strchr(pnt, ',') ) ) pnt++;
			n++;
		}
		if ( 2*(n/2) != n ) {
			sel[n] = molgroup->maxlen - 1;
			n++;
		}
		i = flag = 0;
		for ( j=0; j<molgroup->maxlen; j++ ) {
			if ( j == sel[i] + flag - 1 ) {
				i++;
				flag = 1 - flag;
			}
			molgroup->seqflag[j] *= flag;
		}
	}
*/
	if ( segments.length() )
		seqflag = select_numbers(segments, maxlen);
	
	j = 0;
	for ( i=0; i<maxlen; i++ )
		j += seqflag[i];
	
	if ( verbose & VERB_PROCESS )
		cout << "Alignment positions flagged:    " << j << endl << endl;

	Bresidue_matrix	simat;
	
	simat = get_residue_matrix(simfile);
	
	if ( do_pair ) {
		string		seq1 = seqs[0].sequence();
		string		seq2 = seqs[1].sequence();
		pair<string,string> seq_pair = seq_pair_align(seq1, seq2, gapopen, gapextend, simat);
		seqs[0].sequence(seq_pair.first);
		seqs[1].sequence(seq_pair.second);
	}
		
	if ( consolidate_gaps ) sequence_consolidate_gaps(seqs);
	
	Matrix		mat;
    if ( do_id ) mat = sequence_aligned_identity(seqs, seqflag);

	if ( do_sim ) mat = sequence_aligned_similarity(seqs, seqflag, sim_threshold, simat);

	if ( len2 ) sequence_select(seqs, len1, len2);

	if ( ref ) sequence_select(seqs, mat, ref, cutoff);
	
	if ( delseq ) sequence_delete(seqs, mat);

	if ( do_profile ) profile = sequence_aligned_profile(seqs);
    
	if ( do_hp ) sequence_aligned_hydrophobicity(seqs, seqflag, window, rep_threshold, propfile, psfile);
    
	if ( do_info ) sequence_aligned_information(seqs, seqflag, window, psfile);
	
	if ( optind < argc )
		write_sequence(argv[optind], seqs);

	if ( matfile.length() && mat.rows() )
		mat.write(matfile);

	if ( imgfile.length() && mat.rows() ) {
		Bimage*		pimg = new Bimage(mat, 1);
		pimg->statistics();
		write_img(imgfile, pimg, 0);
		delete pimg;
	}
         
	timer_report(ti);
	
	bexit(0);
}

