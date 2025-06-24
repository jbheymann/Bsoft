/**
@file	seq_util.cpp
@brief	Sequence utility functions 
@author 	Bernard Heymann
@date	Created: 20001029
@date	Modified: 20250601
**/
 
#include "rwgencode.h"
#include "seq_analysis.h"
#include "seq_util.h" 
#include "linked_list.h"
#include "utilities.h" 

#include <map>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
  
/* One-letter to 3-letter mapping */
const map<char, string> res_code = {
        {'-', "GAP"}, 
        {'*', "UNK"}, 
        {'A', "ALA"}, 
        {'B', "ASX"}, 
        {'C', "CYS"}, 
        {'D', "ASP"}, 
        {'E', "GLU"}, 
        {'F', "PHE"}, 
        {'G', "GLY"}, 
        {'H', "HIS"}, 
        {'I', "ILE"}, 
        {'K', "LYS"}, 
        {'L', "LEU"}, 
        {'M', "MET"}, 
        {'N', "ASN"}, 
        {'P', "PRO"}, 
        {'Q', "GLN"}, 
        {'R', "ARG"}, 
        {'S', "SER"}, 
        {'T', "THR"}, 
        {'V', "VAL"}, 
        {'W', "TRP"}, 
        {'Y', "TYR"}, 
        {'Z', "GLX"}, 
        {'X', "UNK"}
} ;

/**
@brief 	Returns the map of one to three character residue codes.
@return map<char, string>		lookup map.
**/
map<char, string>	get_res_code_1_3()
{
	return res_code;
}

/**
@brief 	Returns the map of three to one character residue codes.
@return map<string, char>		lookup map.
**/
map<string, char>	get_res_code_3_1()
{
	map<string, char>	m;
	
	for ( auto c: res_code )
		m[c.second] = c.first;

	return m;
}

/**
@brief 	Get the maxmimum sequences length.
@param 	seqs			set of sequences.
@return long				maximum length.
**/
long		sequence_maximum_length(vector<Bsequence>& seqs)
{
	long		maxlen(0);
	
	for ( auto seq: seqs )
		if ( maxlen < seq.length() ) maxlen = seq.length();

	return maxlen;
}

/**
@brief 	Removes redundant gaps from an alignment.
@param	seqs			set of sequences.
@return int 				0.

	All positions in an alignment with only gaps are removed.

**/
int			sequence_consolidate_gaps(vector<Bsequence>& seqs)
{
	long		nmol(seqs.size());
	long		i, n, m;
	
	long		maxlen = sequence_maximum_length(seqs);
	vector<int>	seqflag(maxlen, 1);
	
	for ( m=i=0; i<maxlen; i++ ) {
		n = 0;
		for ( auto seq: seqs ) {
			if ( i >= seq.length() ) n++;
			else if ( seq.sequence()[i] == '-' ) n++;
		}
		if ( n >= nmol ) {
			seqflag[i] = 0;
			m++;
		}
	}
	
	m = maxlen - m;

	string		newseq(m, '\0');
	
	for ( auto& seq: seqs ) {
		for ( n=i=0; i<maxlen; i++ ) {
			if ( seqflag[i] )
				newseq[n++] = seq.sequence()[i];
		}
		seq.sequence(newseq);
	}
	
	return 0;
}

/**
@brief 	Shows all molecular sequences.
@param 	&seqs 		set of sequences.
@return int				0
**/
int			sequence_show(vector<Bsequence>& seqs)
{
//	long			i, maxnum;
//	Bstring			seq;
	
	if ( verbose & VERB_LABEL )
		cout << endl << "Sequences:" << endl << endl;
	
	for ( auto seq: seqs ) {
		cout << seq.identifier() << tab << seq.length() << endl;
		cout << seq.sequence() << endl;
	}
	
/*	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		cout << "Molecule " << mol->id << ": " << mol->nbase << " nucleotides, "
			<< mol->nres << " residues" << endl;
		if ( mol->naseq.length() ) cout << mol->naseq << endl;
		if ( mol->res ) {
			for ( maxnum=mol->nres, res = mol->res; res; res = res->next )
				if ( maxnum < res->num ) maxnum = res->num;
			seq = Bstring('-', maxnum);
			for ( res = mol->res; res; res = res->next )
				seq[res->num-1] = getcode1(res->type);
			cout << seq << endl;
			seq = 0;
		} else if ( mol->seq.length() ) cout << mol->seq << endl;
		if ( mol->sec ) {
			for ( maxnum=0, sec=mol->sec; sec; sec=sec->next )
				if ( maxnum < sec->last->num ) maxnum = sec->last->num;
			seq = Bstring('-', maxnum);
			for ( sec=mol->sec; sec; sec=sec->next ) {
				if ( sec->type < Strand )
					for ( i=sec->first->num-1; i<sec->last->num; i++ ) seq[i] = 'H';
				if ( sec->type == Strand )
					for ( i=sec->first->num-1; i<sec->last->num; i++ ) seq[i] = 'E';
			}
			cout << seq << endl;
			seq = 0;
		}
		cout << endl;
	}
*/	
	return 0;
}
 
/**
@brief 	Shows the masses of all molecular sequences.
@param 	&seqs 		set of sequences.
@return int				0
**/
int			sequence_mass(vector<Bsequence>& seqs)
{
	long			i, n;
	double			mass;
	string			temp;
	map<char,Bresidue_type>	rtv = get_residue_properties_code1(temp);
		
	if ( verbose & VERB_LABEL )
		cout << "Calculating molecular weights:" << endl << endl;
	
	cout << "Molecule\tMass" << endl;
	for ( auto seq: seqs ) {
		string		s = seq.sequence();
		mass = 0;
		n = seq.length();
		for ( i=0; i<n; i++ )
//			for ( auto rt: rtv ) if ( rt.c == s[i] )
//				mass += rt.mass;
			mass += rtv[s[i]].mass();
		mass -= 18*(n-1);
		cout << seq.identifier() << tab << mass << endl;
	}

	cout << endl;

	return 0;
}

/**
@brief 	Shows the elemental composition of all molecular sequences.
@param 	&seqs 		set of sequences.
@param 	&paramfile		file of residue parameters.
@return vector<double>	array of element numbers: HCNOS
**/
vector<double>	sequence_elements(vector<Bsequence>& seqs, string& paramfile)
{
	long			i, j, n, nat, natr;
	map<char,Bresidue_type>	rtv = get_residue_properties_code1(paramfile);
	vector<double>	comp, el(5,0), elr(5,0), elall(5,0), el1(5,0);
		
	if ( verbose & VERB_LABEL )
		cout << "Calculating the elemental composition:" << endl << endl;
	
	cout << "Molecule\tH\tC\tN\tO\tS\tTotal" << endl;
	for ( auto seq: seqs ) {
		string		s = seq.sequence();
		for ( j=0; j<5; ++j ) el[j] = elr[j] = 0;
		n = s.length();
		for ( i=0; i<n; ++i ) {
//			for ( auto rt: rtv ) if ( rt.c == s[i] )
//				for ( j=0; j<5; ++j ) el[j] += rt.comp[j];
			comp = rtv[s[i]].composition();
			for ( j=0; j<5; ++j ) el[j] += comp[j];
		}
		cout << seq.identifier();
		for ( j=nat=natr=0; j<5; ++j ) {
			nat += el[j];
			natr += elr[j];
			elall[j] += el[j];
			cout << tab << el[j];
		}
		cout << tab << nat << endl;
	}

	cout << "Total";
	for ( j=n=0; j<5; n+=elall[j], ++j ) cout << tab << elall[j];
	cout << tab << n << endl << "Percentage:";
	for ( j=0; j<5; ++j ) cout << tab << elall[j]*100.0/n;
	cout << endl << endl;

	return elall;
}
 

/**
@brief 	Complements all nucleotide sequences.
@param 	&seqs 		set of sequences.
@return int 			0.

	Search through a list of 1-3 mappings for the desired letter.

**/
int 		sequence_complement_all(vector<Bsequence>& seqs)
{
	if ( verbose )
		cout << endl << "Complementing all nucleotide sequences" << endl << endl;
	
	for ( auto& seq: seqs ) 	// Complement all sequences
		complement_sequence(seq.sequence());
	
	return 0;
}

/**
@brief 	Translates all nucleotide sequences to protein sequences.
@param 	&seqs  set of sequences.
@param 	frame			the frame for translation.
@param 	&gcfile			file with genetic code.
@return vector<Bsequence>	translated sequences.

	Each nucleic acid sequence in set of sequences is translated to the
	protein sequence.

**/
vector<Bsequence>	sequence_translate_all(vector<Bsequence>& seqs, int frame, string& gcfile)
{
	vector<Bsequence>	nuseqs;
	map<string,char>	gencode = get_genetic_code(gcfile);
	
	if ( verbose & VERB_LABEL )
		cout << endl << "Translating all nucleotide sequences to amino acid sequences in frame " << frame << endl << endl;
	
	for ( auto seq: seqs ) {
		string		ns = sequence_translate(seq.sequence(), frame, gencode);
		Bsequence	nuseq(seq.identifier());
		nuseq.sequence(ns);
		nuseqs.push_back(nuseq);
	}
	
	return nuseqs;
}

/**
@brief 	Finds a nucleotide sequence.
@param 	&seqs  		set of sequences.
@param 	&seq		sequence to find.
@return long 			position.
**/
/*long 		seq_find_dna(vector<Bsequence>& seqs, Bstring& seq)
{
	Bmolecule*		mol;
	Bstring			sseq = seq.upper();
	Bstring			mol_naseq;
	long			i, offset(0);
	
	if ( verbose )
		cout << "Searching for: " << sseq << endl << endl;
	
	for ( mol = molgroup->mol; mol; mol = mol->next ) {
		mol_naseq = mol->naseq.upper();
		if ( mol_naseq.length() ) {
			offset = mol_naseq.find(sseq, 0);
			if ( offset > -1 ) {
				cout << mol->id << ": " << offset << endl;
				for ( i=0; i<offset; i++ ) cout << " ";
				cout << sseq << endl << mol_naseq << endl;
			} else {
				complement_sequence(mol_naseq);
				offset = mol_naseq.find(sseq, 0);
				if ( offset > -1 ) {
					cout << mol->id << ": " << offset << " (complement)" << endl;
					for ( i=0; i<offset; i++ ) cout << " ";
					cout << sseq << endl << mol_naseq << endl;
				} else {
					cout << mol->id << ": not found" << endl;
				}
			}
		} else
			cout << "Molecule " << mol->id << ": Nucleic acid sequence not found" << endl;
	}
	
	return offset;
}*/

/**
@brief 	Finds a query sequence in a set of sequences.
@param 	&seqs 		set of sequences.
@param 	&findseq	sequence to find.
@return long 			position.
**/
long 		sequence_find(vector<Bsequence>& seqs, string& findseq)
{
	string			sseq = to_upper(findseq);
	long			offset(0);
	
	if ( verbose )
		cout << "Searching for: " << sseq << endl << endl;
	
	for ( auto seq: seqs ) {
		string		s = to_upper(seq.sequence());
		if ( s.length() ) {
			offset = s.find(sseq, 0);
			if ( offset > -1 ) cout << seq.identifier() << ": " << offset << endl;
			else cout << seq.identifier() << ": not found" << endl;
		}
	}
	
	return offset;
}

/**
@brief 	Finds the coding region for an amino acid sequence.
@param 	&seqs 			set of sequences.
@param 	&findseq		sequence to find.
@param 	seqlenmin		sequence length minimum.
@param 	seqlenmax		sequence length maximum.
@param 	side1			preceding sequence length to include.
@param 	side2			succeeding sequence length to include.
@param 	threshold		threshold for reporting possible hits.
@param 	&gcfile			file with genetic code.
@return string 			coding sequence.

	All molecules in the group are searched in all 6 possible frames.

**/
string		seq_find_protein_in_dna(vector<Bsequence>& seqs, string& findseq, int seqlenmin, int seqlenmax, 
				int side1, int side2, double threshold, string& gcfile)
{
	long			n = findseq.length();
	if ( n < 1 ) return "";
	
	if ( seqlenmin < n/2 ) seqlenmin = n/2;
	if ( seqlenmax <= n ) seqlenmin = 2*n;
	
	int				linelength = 80;
	int 			f, i, j, nres = 0, start, length, frame, term, seqlen = 0;
	int 			bestscore = 0, bestframe = 0;
	int 			bestindex = 0, bestseqlen = 0;
	string			bestfrag;
	Bsequence		bestseq;
	char			format[100];
	snprintf(format, 100, "%%.%ds%c", linelength, '\n');
	
	map<string,char>	gencode = get_genetic_code(gcfile);
	
	if ( verbose ) {
		cout << endl << "Finding an amino acid sequence in translated nucleotide sequences:" << endl;
		cout << "Sequence length range:          " << seqlenmin << " - " << seqlenmax << endl;
		cout << "Score threshold (residues):     " << threshold << endl;
	}
	
	for ( auto seq: seqs ) {
		for ( f=0; f<6; f++ ) { 			// Search all translational frames
			string			nucseq = seq.sequence();
			frame = f;
			if ( f > 2 ) {
				frame = f - 3;
				complement_sequence(nucseq);
			}
			string		aaseq = sequence_translate(nucseq, frame, gencode);
			nres = aaseq.length();
			
			vector<int>	score(nres,0);
			
			for ( i=0; i<=nres - n; i++ ) {	// Loop over starting positions
				if ( i == 0 || aaseq[i] == '*' ) {
					seqlen = 0;
					while ( ( i+seqlen+1 < nres ) && ( aaseq[i+seqlen+1] != '*' ) )
						seqlen++;
				}
				if ( seqlen >= seqlenmin && seqlen <= seqlenmax ) {
					term = 0;
					for ( j=0; j<n && j<nres && !term; j++ ) {
						if ( aaseq[i+j] == findseq[j] ) score[i]++;
						if ( aaseq[i+j] == '*' ) term = 1;
					}
					if ( !term && threshold > 0 && score[i] >= threshold ) {
						if ( verbose )
							cout << "Frame: " << f << "  Position: " << 3*i + frame << 
								"   Length: " << seqlen << "   Score: " << score[i] << 
								" (" << score[i]*100.0/findseq.length() << " %)" << endl;
						for ( j=0; j<n; j+= linelength )
							printf(format, &aaseq[i + j]);
					}
					if ( !term && bestscore < score[i] ) {
						bestscore = score[i];
						bestseqlen = seqlen;
						bestframe = f;
						bestindex = 3*i + frame;
						start = (i - side1)*3 + frame;
						while ( start < 0 ) start += 3;
						length = (n + side1 + side2)*3;
						if ( length + start > nucseq.length() )
							length = nucseq.length() - start;
						bestfrag = nucseq.substr(start, length);
						bestseq = seq;
					}
				}
			}
		}
	}
	
	cout << "Best fragment =" << bestfrag << endl;
	
	if ( bestseq.length() && verbose ) {
		cout << endl << "Best sequence score:            " << bestscore << 
			" (" << bestscore*100.0/findseq.length() << "%)" << endl;
		cout << "Sequence:                       " << bestseq.identifier() << endl;
		cout << "Length:                         " << bestseqlen << endl;
		if ( bestframe < 3 ) cout << "Frame:                          " << bestframe << endl;
		else cout << "Frame:                          " << bestframe - 3 << " (complement)" << endl;
		cout << "Index:                          " << bestindex << endl;
		cout << "Nucleotide sequence:" << endl;
		for ( i=0; i<bestfrag.length(); i+= linelength )
			printf(format, bestfrag.c_str() + i);
		string		aaseq = sequence_translate(bestseq.sequence(), 0, gencode);
		cout << "Amino acid sequence:" << endl;
		for ( i=0; i<aaseq.length(); i+= linelength )
			printf(format, aaseq.c_str() + i);
		cout << "Requested sequence:" << endl;
		for ( i=0; i<findseq.length(); i+= linelength )
			printf(format, findseq.c_str() + i);
	}
	cout << endl;
	
	return bestseq.sequence();
}

/**
@brief 	Converts a one-letter amino acid designation to the three-letter equivalent.
@param 	c			the desired amino acid code letter
@param 	*cod		the corresponding three-letter code
@return int			0.

	Search through a list of 1-3 mappings for the desired letter.

**/
int			getcode3(char c, char* cod) 
{ 
	c = toupper(c); 
	
	string		rc(res_code.at(c));
	for ( int i=0; i<3; ++i ) cod[i] = rc[i];
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG getcode3: c=" << c << " cod=" << cod << endl;
		 
	return 0; 
} 

string		getcode3(char c) 
{
	return res_code.at(c);
}

 
/**
@brief 	Converts a three-letter amino acid designation to the one-letter equivalent.
@param 	*acode 		the desired amino acid three-letter code
@return char		the corresponding one-letter code

	Search through a list of 3-1 mappings for the desired three-letter code.

**/
char		getcode1(char* acode) 
{ 
	size_t			i;
	char			c(' ');
	 
	for ( i=0; i<3 && i<strlen(acode); i++ ) acode[i] = toupper(acode[i]);
	acode[i] = 0;
	
	for ( auto it=res_code.begin(); it!=res_code.end(); ++it )
		if ( it->second == acode ) c = it->first;
	 
	return c;
} 
 
/**
@brief 	Complements a nucleotide sequence in place.
@param 	&nucseq		nucleotide sequence to be translated.
@return int 		0.
**/
int 		complement_sequence(string& nucseq)
{
	int 		i, j;
	char		nuc;
	
	j = nucseq.length();
	for ( i=0; i<(j+1)/2; i++ ) {
		nuc = get_complement(nucseq[j-i-1]);
		nucseq[j-i-1] = get_complement(nucseq[i]);
		nucseq[i] = nuc;
	}
	
	return 0;
}

/**
@brief 	Get the Watson-Crick complement of a nucleotide base.
@param 	nuc			nucleotide.
@return char 		complementing nucleotide.
**/
char		get_complement(char nuc)
{
	char		nunuc = 'X';
	
	switch ( nuc ) {
		case 'A': nunuc = 'T'; break;
		case 'G': nunuc = 'C'; break;
		case 'C': nunuc = 'G'; break;
		case 'T': nunuc = 'A'; break;
		case 'U': nunuc = 'A'; break;
		default: break;
	}

	return nunuc;
}

/**
@brief 	Translates a nucleotide sequence to a protein sequence.
@param 	&nucseq		nucleotide sequence to be translated.
@param 	frame		coding frame.
@param 	&gencode	genetic code: array of amino acids.
@return string 		translated protein sequence.
**/
string 		sequence_translate(string& nucseq, long frame, map<string,char>& gencode)
{
	long 			i, j;
	long 			nnuc = nucseq.length();
	long 			nres = (nnuc - frame)/3;
	string			aaseq(nres, ' '), codon;
	
	for ( i=0, j=frame; i<nres && j<nnuc; i++, j+=3 )
		aaseq[i] = gencode[nucseq.substr(j,3)];
	
	if ( i < nres ) aaseq[i] = 0;
	
	return aaseq;
}


