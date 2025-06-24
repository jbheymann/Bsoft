/**
@file	rwGenBank.cpp
@brief	Library routines to read and write GenBank sequence files
@author Bernard Heymann
@date	Created: 20030308
@date	Modified: 20250509
**/

#include "rwGenBank.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads GenBank format sequence files.
@param	*molgroup molecule group.
@return vector<Bsequence>	set of sequences.

	GenBank format:
	LOCUS       NP_039940                284 aa            linear   VRL 13-FEB-2003
	DEFINITION  UL6 [Human herpesvirus 5].
	ACCESSION   NP_039940
	ORIGIN      
	        1 mhakmngwag vrlvthclnt rsrtyvalnm lafartprgv psclfnkvwv sryalvlilm
	       61 vcasesstsw avtsnrlpnc stitttagqd aelhgpapls cnvtqwgrye ngstpvlwct
		  121 lwgsrtrvsl ghrvafgcsw ktffiynvse ssggty
	//

**/
vector<Bsequence>	readGenBank(string& filename)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG readGenBank: filename=" << filename << endl;
	
	vector<Bsequence>	seqs;

    // Open seq file read only
	ifstream		fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return seqs;
	}
    
    string			s, seqstr, seqname;
 	long			i, nid(0), nseq(0), m(0);
	int				seq_flag(0);
	Bsequence		seq;
	
//	int 		gaps = 1;								// Retain gaps
//	if ( strstr(molgroup->select, "NOGAP") ) gaps = 0;	// Omit gaps
	
	// Pass to read the data
	while ( !fseq.eof() ) {
		getline(fseq, s);
		if ( s.length() ) {
			vector<string>	vs = split(s);
			if ( vs[0] == "LOCUS" ) {
				if ( s.find("bp") != string::npos || s.find("DNA") != string::npos )
					seq.type("DNA");
				else 
					seq.type("Protein");
				seqname = vs[1];
				seq_flag = 1;
				nid++;
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG readGenBank: type=" << seq.type() << endl;
			} else if ( vs[0] == "ORIGIN" ) {
				seq_flag++;
			} else if ( vs[0] == "//" ) {
				seqstr = seqstr.substr(0,m);
				seq.identifier(seqname);
				seq.sequence(seqstr);
				seqs.push_back(seq);
				seq_flag = 0;
				seqstr.clear();
				m = 0;
				nseq++;
			} else if ( seq_flag ) {
				s = s.substr(10);
				seqstr += s;
	            for ( i=0; i<s.length(); i++ ) {
    	            if ( isalpha(s[i]) ) {
        	            seqstr[m] = toupper(s[i]);
                	    m++;
	                }
					if ( ( s[i] == '-' ) || ( s[i] == '.' ) ) {
						seqstr[m] = '-';
						m++;
					}
                }
			}
		}
    }
    
	fseq.close();
	
	if ( nseq != nid ) {
		cerr << "Number of identifiers (" << nid << ") don't agree with number of sequences (" << nseq << ")" << endl;
		return seqs;
	}
	
    return seqs;
}


/**
@brief 	Writes GenBank format sequence files.
@param 	&filename	sequence file name.
@param	seqs	 		set of sequences.
@return int 				number of sequences written (<0 if writing failed).
**/
int 	writeGenBank(string& filename, vector<Bsequence> seqs)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG writeGenBank: filename=" << filename << endl;
	
    ofstream        fseq(filename.c_str());
	if ( fseq.fail() ) {	
		cerr << "Error: File " << filename << " not opened!" << endl;
		return -1;
	}

    int				nseq(0);
    long			j, k, n;
    
 	for ( auto seq: seqs ) {
		n = seq.length();
        fseq << "LOCUS       " << seq.identifier() << endl;
//		fseq << mol->id << "," << mol->nres << " bases, 0 checksum" << endl;
        fseq << "ORIGIN" << endl;
        for ( j=0; j<n; j+=60 ) {
			fseq << setw(9) << right << j+1;
			for ( k=j; k<j+60 && k<n; k+=10 )
				fseq << " " << seq.sequence().substr(k,10);
            fseq << endl;
        }
        fseq << "//" << endl << endl;
        nseq++;
    }
    
	fseq.close();
    
    return nseq;
}

