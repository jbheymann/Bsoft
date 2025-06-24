/**
@file	seq_analysis.h 
@brief	Header file for sequence analysis functions
@author 	Bernard Heymann 
@date	Created: 19990123
@date	Modified: 20250510
**/

#include "Bsequence.h"
#include "rwresprop.h"
#include "Complex.h"
#include "Matrix.h"
#include "utilities.h"
  
// Function prototypes 
vector<int>	sequence_limit(vector<Bsequence>& seqs, string& refseq);
Matrix	 	sequence_aligned_identity(vector<Bsequence>& seqs, vector<int> seqflag);
Matrix	 	sequence_aligned_similarity(vector<Bsequence>& seqs, vector<int> seqflag, double threshold, Bresidue_matrix& simat);
long		sequence_select(vector<Bsequence>& seqs, long minlen, long maxlen);
long		sequence_select(vector<Bsequence>& seqs, Matrix mat, long ref, double cutoff);
long		sequence_delete(vector<Bsequence>& seqs, Matrix mat);
string		sequence_aligned_profile(vector<Bsequence>& seqs);
int 	 	sequence_aligned_information(vector<Bsequence>& seqs, vector<int> seqflag, int window, string& psfile);
int 	 	sequence_aligned_hydrophobicity(vector<Bsequence>& seqs, vector<int> seqflag,
				int window, double threshold, string& hphobfile, string& psfile);
vector<Complex<float>>	sequence_frequency_analysis(long win, long start, long end, vector<double>& data);
Matrix		sequence_correlated_mutation(vector<Bsequence>& seqs, vector<int> seqflag,
					string& refseqid, double cutoff, Bresidue_matrix& simat);
