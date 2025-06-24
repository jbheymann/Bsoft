/**
@file	seq_util.h 
@brief	Header file for sequence utilities 
@author 	Bernard Heymann 
@date	Created: 20001029
@date	Modified: 20250601
**/
 
#include "Bsequence.h"

// Function prototypes 
map<char, string>	get_res_code_1_3();
map<string, char>	get_res_code_3_1();
long		sequence_maximum_length(vector<Bsequence>& seqs);
int			sequence_consolidate_gaps(vector<Bsequence>& seqs);
//int			seq_from_residues(vector<Bsequence>& seqs);
int			sequence_show(vector<Bsequence>& seqs);
int			sequence_mass(vector<Bsequence>& seqs);
vector<double>	sequence_elements(vector<Bsequence>& seqs, string& paramfile);
int 		sequence_complement_all(vector<Bsequence>& seqs);
vector<Bsequence>	sequence_translate_all(vector<Bsequence>& seqs, int frame, string& gcname);
long 		sequence_find_dna(vector<Bsequence>& seqs, string& seq);
long 		sequence_find_protein(vector<Bsequence>& seqs, string& seq);
string		sequence_find_protein_in_dna(vector<Bsequence>& seqs, string& seq, int seqlenmin, int seqlenmax, 
				int side1, int side2, double threshold, string& gcfile);
int			getcode3(char c, char* cod);
string		getcode3(char c);
char		getcode1(char* acode); 
int 		complement_sequence(string& nucseq);
char		get_complement(char nuc);
string 		sequence_translate(string& nucseq, long frame,  map<string,char>& gencode);

