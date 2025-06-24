/**
@file	seq_align.h
@brief	Header file for functions to generate and analyze dot plots
@author 	Bernard Heymann
@date	Created: 20001029
@date	Modified: 20250509
**/

#include "rwresprop.h"

// Function prototypes 
pair<string,string>	seq_pair_align(string seq1, string seq2, double gapopen, 
				double gapextend, Bresidue_matrix& simat);
int			seq_find_best_offset(string seq1, string seq2, long& nres, Bresidue_matrix& simat);
Matrix		seq_dot_plot(string seq1, string seq2, Bresidue_matrix& simat);
Matrix		seq_dot_plot_mov_avg(Matrix dot_plot, int window);
int			seq_dot_plot_interpret(Matrix dot_plot);
double		seq_dot_plot_best_segments(Matrix dot_plot);

