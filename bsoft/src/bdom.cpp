/**
@file	bdom.c
@brief	A program to generate protein domain diagrams
@author Bernard Heymann
@date	Created: 20110412 
@date	Modified: 20250602
**/

#include "ps_plot.h" 
#include "ps_sequence.h" 
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

vector<Bgroup>	domain_read(string filename);

// Usage assistance
const char* use[] = {
" ",
"Usage: bdom [options] input.txt output.ps",
"-----------------------------------------",
"Generates protein domain diagrams.",
" ",
"Actions:",
"-numbers                 Show residue numbers.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-name vp25               Molecule name.",
"-width 1374              Width of diagram (in sequence numbers).",
"-height 60               Height of diagram.",
" ",
//"Input:",
//"-properties file.star    Property file (default = internal values).",
//" ",
//"Output:",
//"-Image file.img          Image output filename.",
//" ",
NULL
};

int 		main(int argc, char **argv)
{
    // Initialize variables
	bool			numbers(0);			// Flag to show residue numbers
	string			name;				// Molecule name
	long			width(0);			// Width of diagram
	long			height(20);			// Height of diagram
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "numbers" ) numbers = 1;
		if ( curropt->tag == "name" )
			name = curropt->value.str();
		if ( curropt->tag == "width" )
			if ( ( width = curropt->integer() ) < 1 )
				cerr << "-width: A number of residues must be specified!" << endl;
		if ( curropt->tag == "height" )
			if ( ( height = curropt->integer() ) < 1 )
				cerr << "-height: A value must be specified!" << endl;
    }
	option_kill(option);
    
	double		ti;
	if ( verbose & VERB_TIME )
		ti = timer_start();

    // Read the domain specification file
    vector<Bgroup>	doms = domain_read(argv[optind++]);
	if ( doms.size() < 1 ) {
		cerr << "Error: No molecule domains read!" << endl;
		bexit(-1);
	}
	
	if ( name.length() < 1 ) name = "Domains";
	
	string		filename(argv[optind]);
	ps_plot_domains(filename, name, width, height, doms, numbers);
	
	if ( verbose & VERB_TIME )
		timer_report(ti);
	
	bexit(0);
}

vector<Bgroup>	domain_read(string filename)
{
	string			s;
 	long			n(0), start, end;
	RGBA<float> 	color;
	vector<Bgroup>	doms;

   // Open domain file read only
	ifstream		fdom(filename.c_str());
	if ( fdom.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return vector<Bgroup>();
	}
    	
	// Read the domain data
	while ( !fdom.eof() ) {
		getline(fdom, s);
		if ( s.length() ) {
			istringstream	ss(s);
			ss >> start >> end >> color[0] >> color[1] >> color[2] >> color[3];
			s = to_string(++n);
			doms.push_back(Bgroup(s, start, end, color));
		}
	}
	
	fdom.close();
	
	return doms;
}

