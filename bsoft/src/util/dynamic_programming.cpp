/**
@file	dynamic_programming.cpp
@brief	Library functions for dynamic programming.
@author 	Bernard Heymann
@date	Created: 20050622
@date	Modified: 20250509
**/

#include "dynamic_programming.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the scoring matrix in dynamic programming.
@param 	mat			matrix.
@param 	gapopen		gap opening penalty.
@param 	gapextend	gap extension penalty.
@return long 			index of the maximum in the matrix.

	Implementation of the Needleman-Wunsch algorithm.
	The input matrix can be any form of similarity or coincidence matrix, and
	is modified by dynamic programming in preparation for backtracing.

**/
long		dp_matrix_scoring(Matrix mat, double gapopen, double gapextend)
{
	long		m = mat.columns(), n = mat.rows();
	long		t, pt(0);
	long		i, j, k, kmax(0);
	double		value, max, maxall(0);
	
//	if ( verbose & VERB_DEBUG )
		cout << "DEBUG dp_matrix_scoring: rows=" << n << " columns=" << m << endl;
	
	for ( k=j=0; j<n; j++ ) {
		for ( i=0; i<m; i++, k++ ) {
			if ( i==0 || j==0 ) max = mat[j][i];
			else max = mat[j][i] + mat[j-1][i-1];
			value = 0;
			t = 0;
			if ( i > 0 ) {
				if ( pt == 1 ) value = mat[j][i-1] - gapextend;
				else value = mat[j][i-1] - gapopen;
			}
			if ( max < value ) {
				max = value;
				t = 1;
			}
			if ( j > 0 ) {
				if ( pt == 2 ) value = mat[j-1][i] - gapextend;
				else value = mat[j-1][i] - gapopen;
			}
			if ( max < value ) {
				max = value;
				t = 2;
			}
			mat[j][i] = max;
			if ( maxall < max ) {
				maxall = max;
				kmax = k;
			}
			pt = t;
		}
	}
	
//	if ( verbose & VERB_DEBUG )
		cout << "DEBUG dp_matrix_scoring: kmax=" << kmax << endl;

	return kmax;
}

/**
@brief 	Backtraces the scoring matrix in dynamic programming.
@param 	mat			matrix.
@param 	gapopen		gap opening penalty.
@param 	gapextend	gap extension penalty.
@return int* 			alignment array with indices.

	Implementation of the Needleman-Wunsch algorithm.
	The scoring matrix for dynamic programming is backtraced into an integer 
	array holding the indices from the corresponding sequences, or -1 for gaps. 
	The array is double the length of the alignment and holds two sets of 
	indices, the first in the first half and the second in last half.
	The length of the alignment is returned in the pointer of the last argument.

**/
pair<vector<int>, vector<int>>	dp_matrix_backtrace(Matrix mat, double gapopen, double gapextend)
{
	long		m = mat.columns(), n = mat.rows();
	long		t, pt(0);
	long		i, j, h, hbeg(0), hend, imax(0), jmax(0);
	long		tlen = m + n;
	double		max(0), tmax(0), value;
	vector<int>	aln1(tlen, -1);
	vector<int>	aln2(tlen, -1);

//	if ( verbose & VERB_DEBUG )
		cout << "DEBUG dp_matrix_backtrace: rows=" << n << " columns=" << m << endl;
	
	for ( j=0; j<n; j++ ) {
		for ( i=0; i<m; i++ ) {
			if ( tmax < mat[j][i] ) {
				tmax = mat[j][i];
				imax = i;
				jmax = j;
			}
		}
	}
	
	// Fill in the trailing pieces
	hend = m - imax;
	if ( n - jmax > hend ) hend = n - jmax;
	hend = tlen - hend;
	for ( h=hend, i=imax; h<tlen && i<m; h++, i++ ) aln2[h] = i;
	for ( h=hend, j=jmax; h<tlen && j<n; h++, j++ ) aln1[h] = j;
	
	// Trace back to get the central aligned pieces
	for ( h=hend, i=imax, j=jmax; i>0 && j>0; h-- ) {
		mat[j][i] = 2*tmax;
		t = 0;
		max = mat[j-1][i-1];
		if ( pt == 1 ) value = mat[j][i-1] - gapextend;
		else value = mat[j][i-1] - gapopen;
		if ( max < value ) {
			max = value;
			t = 1;
		}
		if ( pt == 2 ) value = mat[j-1][i] - gapextend;
		else value = mat[j-1][i] - gapopen;
		if ( max < value ) {
			max = value;
			t = 2;
		}
		if ( t < 2 ) {
			aln2[h] = i;
			i--;
		}
		if ( t%2 == 0 ) {
			aln1[h] = j;
			j--;
		}
		pt = t;
	}
	mat[j][i] = 2*tmax;
	
	// Fill in the leading pieces
	for ( ; h >= 0; i--, j--, h-- ) {
		if ( i >= 0 ) aln2[h] = i;
		if ( j >= 0 ) aln1[h] = j;
		if ( i >= 0 || j >= 0 ) hbeg = h;
	}
	
	// Generate the return arrays
	long			len = tlen - hbeg;
	vector<int>		align1(len, '-');
	vector<int>		align2(len, '-');
	
	for ( i=0, h=hbeg; i<len; i++, h++ ) {
		align2[i] = aln2[h];
		align1[i] = aln1[h];
	}

	pair<vector<int>, vector<int>>	align_ind = make_pair(align1, align2);
	
	return align_ind;
}
