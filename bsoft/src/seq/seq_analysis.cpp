/**
@file	seq_analysis.cpp
@brief	Analysis of protein sequences 
@author 	Bernard Heymann
@date	Created: 19990123
@date	Modified: 20250510
**/
 
#include "seq_analysis.h"
#include "seq_util.h" 
#include "ps_sequence.h"
#include "moving_average.h"
#include "utilities.h" 
 
// Declaration of global variables
extern int 		verbose;		// Level of output to the screen
extern string 	defaultcode;	// Residue one letter code

double			log_2(double a)
{
	if ( a <= 0 ) return 0;
	else return log(a)/log(2.0L);
}

/**
@brief 	Limits the selection to the reference sequence in an aligned set.
@param 	*molgroup 		the set of sequences.
@param 	&refseq 			reference sequence identifier.
@return long 				number of selected residues.

**/
vector<int>	sequence_limit(vector<Bsequence>& seqs, string& refseq)
{
	long			i;
	
	long			maxlen = sequence_maximum_length(seqs);
	
	vector<int>	seqflag(maxlen, 1);
	
	for ( i=0; i<seqs.size(); ++i )
		if ( seqs[i].identifier().find(refseq) != string::npos ) break;
		
	if ( i >= seqs.size() ) {
		cerr << "Sequence " << refseq << " not found!" << endl;
		return seqflag;
	}
	
	string&		seq = seqs[i].sequence();
	
	for ( i=0; i<seq.length(); ++i ) {
		if ( seq[i] == '-' )
			seqflag[i] = 0;
		else
			seqflag[i] = 1;
	}
	
	return seqflag;
}

/**
@brief 	Calculates the pairwise identities between aligned sequences.
@param 	seqs	 		the set of sequences.
@param 	seqflag	 		boolean array to select sequence positions.
@return Matrix 			the matrix of identities.

The identity between two sequences is defined as:
			   number of identical residues
	identity = ----------------------------
					   overlap
where the overlap is the number of positions with residues in both sequences.

**/
Matrix	 	sequence_aligned_identity(vector<Bsequence>& seqs, vector<int> seqflag)
{
	long   			i, j, k, n(0), nid, overlap, nseq(seqs.size());
	long			idsum(0), idssum(0), overlapsum(0), overlapssum(0);
	double			idavg, overlapavg;

	long			maxlen = sequence_maximum_length(seqs);
	
	if ( seqflag.size() < maxlen )
		seqflag.resize(maxlen, 1);
	
	if ( verbose & VERB_LABEL )
	    cout << "Aligned identity analysis:" << endl;

	if ( verbose & VERB_PROCESS )
		cout << "Seq1\tSeq2\tIdentity\tnID\tOverlap\tName1\tName2" << endl;
	
	// Initialize an image structure to hold the results
	Matrix		mat(nseq, nseq);
	
	for ( i=1; i<nseq; ++i ) {
		string&		seq1 = seqs[i].sequence();
		for ( j=0; j<i; ++j ) {
			string&		seq2 = seqs[j].sequence();
			n++;
			nid = overlap = 0;
			maxlen = seq1.length();
			if ( maxlen > seq2.length() ) maxlen = seq2.length();
			for ( k=0; k<maxlen; k++ ) {
				if ( seqflag[k] > 0 &&
						seq1[k] != '-' && 
						seq2[k] != '-' ) {
					overlap++;
					if ( seq1[k] == seq1[k] ) nid++;
				}
			}
			idsum += nid;
			idssum += nid*nid;
			overlapsum += overlap;
			overlapssum += overlap*overlap;
			if ( overlap ) mat[i][j] = mat[j][i] = nid*1.0/overlap;
			if ( verbose & VERB_PROCESS )
				cout << i+1 << tab << j+1 << tab << fixed << setprecision(3) << setw(8) <<
					mat[i][j] << tab << nid << tab << overlap << tab <<
					seqs[i].identifier() << tab << seqs[j].identifier() << endl;
		}
	}
	
	idavg = idsum*1.0/n;
	overlapavg = overlapsum*1.0/n;
	
	if ( verbose & VERB_LABEL ) {
		cout << "Average identical residues:  " << idavg << " (" << 
			sqrt(idssum*1.0/n - idavg*idavg) << ")" << endl;
		cout << "Average overlap:             " << overlapavg << " (" << 
			sqrt(overlapssum*1.0/n - overlapavg*overlapavg) << ")" << endl << endl;
	}
	
    return mat;
}

/**
@brief 	Calculates the pairwise similarities between aligned sequences.
@param 	seqs	 		the set of sequences.
@param 	seqflag	 		boolean array to select sequence positions.
@param 	threshold 		threshold to accept residues as similar.
@param 	&simat			residue similarity matrix.
@return Matrix			the matrix of similarities.

The similarity between two sequences is defined as:
				   sum(residue similarity)
	similarity   = -----------------------
						  overlap
						  number of residues with similarity > threshold
	fraction similarity = ----------------------------------------------
											overlap
where the overlap is the number of positions with residues in both sequences.
The residue similarity is taken from a residue substitution matrix.
The default substitution matrix is BLOSUM62.

**/
Matrix	 	sequence_aligned_similarity(vector<Bsequence>& seqs, vector<int> seqflag, double threshold, Bresidue_matrix& simat)
{
    long 				i, j, k, thelen, overlap, nseq(seqs.size());
    long				ii, jj;
	double				simsum(0), percentsim(0);
	
	long			maxlen = sequence_maximum_length(seqs);
	
	if ( seqflag.size() < maxlen )
		seqflag.resize(maxlen, 1);
	
	string				code1 = simat.code();
	Matrix				sim_mat = simat.matrix();
    
	if ( verbose & VERB_LABEL ) {
	    cout << "Aligned similarity analysis:" << endl;
	    cout << "Similar residue threshold:      " << threshold << endl;
	    cout << "Number of sequences:            " << nseq << endl << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Seq1\tSeq2\tSimilarity\tFraction\tOverlap\tName1\tName2" << endl;
	
	// Initialize an image structure to hold the results
	Matrix		mat(nseq, nseq);
	
	for ( i=1; i<nseq; ++i ) {
		string&		seq1 = seqs[i].sequence();
		for ( j=0; j<i; ++j ) {
			string&		seq2 = seqs[j].sequence();
            thelen = seq1.length();
            if ( thelen > seq2.length() ) thelen = seq2.length();
            simsum = percentsim = 0.0;
			overlap = 0;
			if ( verbose & VERB_DEBUG )
				cout << "i=" << i << " j=" << j << " thelen=" << thelen << endl;
			for ( k=0; k<thelen; k++ ) {
				if ( seqflag[k] > 0 &&
						( seq1[k] != '-' ) && ( seq2[k] != '-' ) ) {
					ii = code1.find(seq1[k]);
					jj = code1.find(seq2[k]);
					if ( ii > -1 && jj > -1 ) {
	                	simsum += sim_mat[ii][jj];
						if ( sim_mat[ii][jj] >= threshold ) 
							percentsim += 1;
						overlap++;
					}
				}
				if ( verbose & VERB_DEBUG )
					cout << k << tab << simsum << tab << overlap << endl;
			}
			if ( overlap ) {
				mat[i][j] = simsum*1.0/overlap;
				percentsim /= overlap;
			} else mat[i][j] = 0;
			mat[j][i] = mat[i][j];
			if ( verbose & VERB_PROCESS )
				cout << i+1 << tab << j+1 << tab << fixed << setprecision(3) << setw(8) <<
					mat[i][j] << tab << setw(8) << percentsim << tab << overlap << tab <<
					seqs[i].identifier() << tab << seqs[j].identifier() << endl;
        }
	}
	cout << endl;
    
    return mat;
}

/**
@brief 	Selects sequences within a range of lengths.
@param 	seqs			the set of sequences.
@param 	minlen 			minimum length.
@param 	maxlen 			maximum length.
@return long				number of sequences retained.
**/
long		sequence_select(vector<Bsequence>& seqs, long minlen, long maxlen)
{
	long			i, n(0), len;

	if ( verbose )
		cout << "Selecting sequences of length " << minlen << " - " << maxlen << endl;

	for ( auto seq: seqs ) {
		for ( i=len=0; i<seq.length(); ++i )
			if ( seq.sequence()[i] != '-' ) len++;
		if ( len >= minlen && len <=maxlen ) {
			seq.select(1);
			n++;
		} else {
			seq.select(0);
		}
	}

	if ( verbose )
		cout << "Selected:                  " << n << endl << endl;

	return n;
}

/**
@brief 	Selects sequences based on a comparison matrix of aligned sequences.
@param 	seqs			the set of sequences.
@param 	mat 				comparison matrix.
@param 	ref		 		reference sequence number (starting at 1).
@param 	cutoff	 		threshold for selecting sequences.
@return long				number of sequences retained.
**/
long		sequence_select(vector<Bsequence>& seqs, Matrix mat, long ref, double cutoff)
{
	if ( ref > mat.rows() ) return mat.rows();
	
	long			i(ref-1), j, n(0);
	
	if ( verbose ) {
		cout << "Selecting based on sequence     " << ref << ":" << endl;
		cout << "Sequence ID:                    " << seqs[i].identifier() << endl;
		cout << "Cutoff:                         " << cutoff << endl << endl;
	}
	
	for ( j = 0; j < mat.rows(); ++j )
		if ( i != j && mat[i][j] < cutoff ) {
			seqs[j].select(0);
		} else {
			seqs[j].select(1);
			n++;
		}

	if ( verbose )
		cout << "Selected:                  " << n << endl << endl;
	
	return n;
}

/**
@brief 	Deletes non-selected sequences and corresponding elelments of a comparison matrix.
@param 	seqs			the set of sequences.
@param 	mat 				comparison matrix.
@return long				number of sequences retained.
**/
long		sequence_delete(vector<Bsequence>& seqs, Matrix mat)
{
	long			i(0);
	vector<int>	del(seqs.size(), 0);

	for ( auto seq: seqs )
		if ( seq.select() == 0 ) {
			++i;
			del[i] = 1;
		}

	if ( verbose )
		cout << "Deleting " << i << " sequences" << endl << endl;
	
//	cout << mat << endl;

	if ( mat.rows() )
		for ( i = mat.rows() - 1; i >= 0; --i )
			if ( del[i] ) mat = mat.delete_row_column(i);
	
	for ( i = seqs.size() - 1; i >= 0; --i )
		if ( del[i] ) seqs.erase(seqs.begin() + i);

	if ( verbose ) {
		cout << "Number of sequences retained:   " << seqs.size() << endl;
		cout << "Number of rows retained:        " << mat.rows() << endl << endl;
	}
	
	return mat.rows();
}

/**
@brief 	Generates a PROSITE format profile from an aligned set of sequences.
@param 	seqs			the set of sequences.
@return string			profile in PROSITE format.

At each position in the alignment, the number of distinct residue types
are counted. If there are more than 3 residue types represented at a
position, or there is a gap, it is designated as variable by an "x".
The profile finally contains 1-3 residue type possibilities for highly
conserved positions interspersed by variable length gaps.

**/
string		sequence_aligned_profile(vector<Bsequence>& seqs)
{
	bool			notdone;
    long 			i, j, k, x, g;
	char			reslist[3];

	long			maxlen = sequence_maximum_length(seqs);
	
	long			len(3*maxlen);
	string			profseq(len,' ');
	string			profile;
	
	if ( verbose & VERB_LABEL )
	    cout << "Profile:" << endl;
	
	for ( i=k=0; i<maxlen; i++, k+=3 ) {
		notdone = 1;
		for ( j=0; j<3; j++ ) reslist[j] = ' ';
		for ( auto seq: seqs ) {
			if ( seq.sequence()[i] == '-' ) {
				for ( j=k; j<k+3; j++ ) profseq[j] = '-';
				notdone = 0;
			} else {
				for ( j=0; j<3 && reslist[j] != ' ' && seq.sequence()[i] != reslist[j]; j++ ) ;
				if ( j<3 ) reslist[j] = seq.sequence()[i];
				else {
					for ( j=k; j<k+3; j++ ) profseq[j] = 'x';
					notdone = 0;
				}
			}
			if ( notdone ) for ( j=0; j<3; j++ ) profseq[k+j] = reslist[j];
		}
	}
	
	for ( i=g=x=k=0; i<maxlen; i++, k+=3 ) {
		if ( profseq[k] == '-' ) g++;
		else if ( profseq[k] == 'x' ) x++;
		else {
			if ( profile.length() ) profile += '-';
			if ( g + x > 0 ) {
				if ( g == 0 ) {
					profile += "x(" + to_string(x) + ")-";
				} else {
					profile += "x(" + to_string(x) + "," + to_string(g) + ")-";
				}
			}
			if ( profseq[k+1] == ' ' ) {
				profile += profseq[k];
			} else {
				profile += '[';
				for ( j=0; j<3 && profseq[k+j] != ' '; j++ ) profile += profseq[k+j];
				profile += ']';
			}
			g = x = 0;
		}
	}
	
	if ( verbose )
		cout << profile << endl;
	
    return profile;
}

/**
@brief 	Calculates the sequence logo representation for an alignment.
@param 	seqs	 		the set of sequences.
@param 	seqflag	 		boolean array to select sequence positions.
@param 	window			window for calculating the moving average.
@param 	&psfile			the postscript file name.
@return int 				0.

The information content of each position in an alignment is calculated
as:
	information = log_2(n) - sum(pi * log_2(pi) )
					 fi
	pi          =  -------
				   sum(fi)
	fi          =  frequency of residue type i at this position
	n           =  sum(fi) if sum(fi) < 20, otherwise n = 20
A moving average of the information is calculated over a given window
to smooth the resultant data.
The sequence logo representation for the occurrence of every residue type
at every position is generated and written into a postscript file.

**/
int 	 	sequence_aligned_information(vector<Bsequence>& seqs, vector<int> seqflag, int window, string& psfile)
{
	long			i, j, k, ii, n, nseq(seqs.size());
 	double       	fgap;
	
	long			maxlen = sequence_maximum_length(seqs);
	
	if ( seqflag.size() < maxlen )
		seqflag.resize(maxlen, 1);

	long			npos(0);
	for ( i=0; i<maxlen; i++ ) npos += seqflag[i];

	string			code1(defaultcode);
	long			nres(code1.size());
    
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_aligned_information: code=" << code1 << endl;
    
	
	long   			igap = code1.find_first_of('-');
	vector<double>	fsum(npos,0);
	vector<double>	info(npos,0);
    vector<double>	freq(nres*npos);
	string			pattern(nres*npos,' ');
	string			seq(npos,' ');
	
    if ( verbose & VERB_LABEL ) {
		cout << "Information content analysis:" << endl;
		cout << "Sequences:                      " << nseq << endl;
		cout << "Alignment positions:            " << npos << endl;
		cout << "Moving average window:          " << window << endl;
	}
	
    // Calculating the frequency of each residue type occurring at each
    // position in the alignment
    for ( auto seq: seqs ) {
    	string&		s = seq.sequence();
    	j = k = 0;
    	while ( j < s.length() ) {
			if ( seqflag[j] ) {
				ii = code1.find_first_of(s[j]);
    	        if ( ii > -1 ) freq[nres*k+ii] += 1;
				k++;
			}
    	    j++;
    	}
    }
	
    // Calculating information content for a position in the alignment
	multimap<double,char> 	ri;
	double					v;
	k = 0;
    for ( j=0; j<maxlen; j++ ) {
        if ( seqflag[j] ) {
	    	fgap = 0;
	    	ri.clear();
    		for ( i=0; i<nres; i++ ) {
				if ( i == igap )
					fgap += freq[nres*k+i];
				else
    		    	fsum[k] += freq[nres*k+i];
			}
			info[k] = 0;
			if ( fsum[k] > 0 ) {
    			for ( i=0; i<nres; i++ ) if ( i != igap ) {
					v = freq[nres*k+i]*log_2(freq[nres*k+i])/fsum[k];
    	    		ri.insert(pair<double,char>(v,code1[i]));
    	    		info[k] += v;
				}
				if ( fsum[k] > 20 )
					info[k] += log(20.0/fsum[k])/log(2.0);
			}
			i = 0;
			for ( auto it: ri ) {
				pattern[nres*k+i] = it.second;
				freq[nres*k+i] = it.first;
				i++;
			}
			seq[k] = seqs[0].sequence()[j];
			k++;
			if ( verbose & VERB_FULL )
				cout << j << tab << fgap << endl;
		}
    }
	
	// Calculate a moving average
	vector<double>		movavg = moving_average(info, window);
	
	// Print out the list of results
    if ( verbose & VERB_LABEL )
		cout << "  Pos\tBits\tMovAvg\tnSeq\tPattern" << endl;
    for ( k=0; k<npos; k++ ) {
    	cout << k+1 << tab << info[k] << tab << movavg[k] << tab << fsum[k] << tab;
		for ( i=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 )
			cout << pattern[nres*k+i];
		cout << endl;
    }
	
    for ( k=0; k<npos; k++ ) {
		if ( k > 0 ) cout << "-";
		for ( i=n=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 ) n++;
		if ( n < 1 ) {
			cout << "x";
		} else {
			if ( n > 1 ) cout << "[";
			for ( i=0; i<nres; i++ ) if ( freq[nres*k+i] > 0.1 )
				cout << pattern[nres*k+i];
			if ( n > 1 ) cout << "]";
		}
    }
	cout << endl;
	
	vector<Complex<float>>	per = sequence_frequency_analysis(16, 4, 4, info);
	
	vector<Complex<float>>	per_avg = moving_average_complex(per, window);
	
	Bstring				title("Information"), pf(psfile);
	if ( psfile.length() )
		ps_seq_info(pf, title, nres, movavg, fsum, per_avg,
			freq, pattern);
	
	return 0;
}

/**
@brief 	Calculates the average hydrophobicity at every position in an alignment.
@param 	seqs	 		the set of sequences.
@param 	seqflag	 		boolean array to select sequence positions.
@param 	window			moving average window.
@param 	threshold		fraction of sequences with a residue in a position.
@param 	&hphobfile		parameter file.
@param 	&psfile			postscript output file.
@return int 	 			0.

The default hydrophobicity scale is the GES scale.

**/
int 	 	sequence_aligned_hydrophobicity(vector<Bsequence>& seqs, vector<int> seqflag,
				int window, double threshold, string& hphobfile, string& psfile)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_aligned_hydrophobicity: " << hphobfile << endl;

	long			maxlen = sequence_maximum_length(seqs);
	
	if ( seqflag.size() < maxlen )
		seqflag.resize(maxlen, 1);

	long					i, j, nseq(seqs.size());	
	map<char,Bresidue_type>	rtv = get_residue_properties_code1(hphobfile);
//	string					code1;
//	for ( auto rt: rtv ) code1.push_back(rt.c);
	
    if ( verbose & VERB_LABEL )
		cout << "Average aligned hydrophobicity:" << endl;
    
	vector<double>	ncnt(maxlen,0);
	vector<double>	avgHP(maxlen,0);
	vector<int>		HPseg(maxlen,0);
	
	for ( auto seq: seqs ) {
		string&			s = seq.sequence();
        for ( j=0; j<s.length(); j++ ) {
			if ( s[j] != '-' ) {
//				if ( ( i = code1.find(s[j]) ) != string::npos ) {
//					avgHP[j] += -rtv[i].hphob; 	// Note: Negative for scales of hydrophilicity
//					ncnt[j]++;
//				}
				avgHP[j] += -rtv[s[j]].hydrophobicity();
				ncnt[j]++;
			}
        }
    }

    for ( j=0; j<maxlen; j++ )
        if ( ncnt[j] > 0 ) avgHP[j] /= ncnt[j];
    
    if ( verbose & VERB_PROCESS )
		cout << "nRes\tavgHP\tintHP\tnSeq\tRes" << endl;
	
	double		sum(0);
	i = 0;
    for ( j=0; j<maxlen; j++ ) {
        if ( seqflag[j] ) {
			i++;
			sum += avgHP[j];
			cout << i << tab << avgHP[j] << tab << sum << tab << ncnt[j] << tab << seqs[0].sequence()[j] << endl;
		}
	}
	
    if ( verbose & VERB_PROCESS )
		cout << endl << "Segment\tStart\tEnd\tLength\tsumHP" << endl;
	
	int 		k, m, sum_flag(0), start[1000], end[1000], newstart, newend;
	double		max, sumHP[1000];
	
	if ( threshold < 0.5 ) threshold = 0.5;
	if ( threshold > 1 ) threshold = 1;
	threshold *= nseq;
	i = 0;
	sumHP[0] = 0;
    for ( j=1; j<maxlen; j++ ) {
        if ( seqflag[j] ) {
			if ( ncnt[j] >= threshold && ncnt[j-1] < threshold ) {
				if ( sumHP[i] > 20 ) i++;
				start[i] = j;
				sum_flag = 1;
				sumHP[i] = 0;
			}
			if ( sum_flag && ncnt[j] < threshold ) {
				end[i] = j - 1;
				sum_flag = 0;
				newstart = start[i];
				newend = end[i];
				max = sumHP[i];
				for ( k=start[i]; k<end[i]; k++ ) {
					sum = 0;
					for ( m=k; m<=end[i]; m++ ) {
						sum += avgHP[m];
						if ( max < sum ) {
							max = sum;
							newstart = k;
							newend = m;
						}
					}
				}
				start[i] = newstart;
				end[i] = newend;
				sumHP[i] = max;
				if ( max > 20 ) {
					for ( k=start[i]; k<=end[i]; k++ ) HPseg[k] = 1;
					cout << i+1 << tab << start[i]+1 << tab << end[i]+1 << tab << 
						end[i] - start[i] + 1 << tab << sumHP[i] << endl;
				}
			}
			if ( sum_flag ) sumHP[i] += avgHP[j];
		}
	}
	cout << endl;
	
	// Calculate a moving average
	vector<double>			movavg = moving_average(avgHP, window);
	
	vector<Complex<float>>	per = sequence_frequency_analysis(16, 4, 4, avgHP);
	
	vector<Complex<float>>	per_avg = moving_average_complex(per, window);
	
	Bstring				title("Hydrophobicity"), pf(psfile);
	if ( psfile.length() )
		ps_seq_hydrophob(pf, title, movavg, HPseg, ncnt, per_avg);
	
	return 0;
}

/**
@brief 	Fourier transforms a vector for frequency analysis.
@param 	win				window size.
@param 	start			start within window.
@param 	end				end within window.
@param 	*data			sequence.
@return int 	 			0.

A brute force Fourier transform is done.

**/
vector<Complex<float>>	sequence_frequency_analysis(long win, long start, long end, vector<double>& data)
{
	if ( start > end ) swap(start, end);
	if ( start < 0 ) start = 0;
	if ( end > win - 1 ) end = win - 1;
	
    if ( verbose )
		cout << "Frequency analysis:" << endl;
       
	if ( verbose & VERB_RESULT ) {
		cout << "Window:                         " << win << endl;
		cout << "Start - end:                    " << start << " - " << end << endl;
	}
	
    long         	j, k, l, length(data.size());
	long			half_win = win/2;
	double			real, imag, amp, phi, dphi;
	
    vector<Complex<float>>	transform_data(length*(end-start+1));
	
    for ( j=0; j<length-win; j++ ) {
		cout << j+half_win;
		for ( k=start; k<=end; k++ ) {
			real = 0;
			imag = 0;
            for ( l=0; l<win; l++ ) {
				phi = MIN2PI*k*l/win;
                real += data[l+j]*cos(phi);
                imag += data[l+j]*sin(phi);
			}
			amp = sqrt(real*real + imag*imag);
            phi = atan2(imag, real);
			if ( amp < 0.001 ) phi = 0;
            dphi = phi - TWOPI*j*k/win;
			dphi = angle_set_negPI_to_PI(dphi);
			transform_data[length*(k-start)+j+half_win] = Complex<float>(amp*cos(dphi), amp*sin(dphi));
			if ( verbose & VERB_RESULT )
				cout << tab << amp << tab << dphi*180/M_PI;
		}
		if ( verbose & VERB_RESULT )
			cout << endl;
    }
	cout << endl;
	
    return transform_data;
}

vector<double>		sequence_aligned_weight(vector<Bsequence>& seqs)
{
	char			m1, m2, igap = '-';
	long   			i, k, l, nseq(seqs.size()), overlap, identity;

	long			maxlen = sequence_maximum_length(seqs);
	
	vector<double>	weight(nseq*nseq, 0);
	
	// Weighting is based on the difference between identity and overlap
	for ( k=1; k<nseq; ++k ) {
		string&		seq1 = seqs[k].sequence();
		for ( l=0; l<k; ++l ) {
			string&		seq2 = seqs[l].sequence();
			overlap = identity = 0;
			for ( i=0; i<maxlen; i++ ) {
				m1 = seq1[i];
				m2 = seq2[i];
				if ( m1 != igap && m2 != igap ) {
					overlap++;
					if ( m1 == m2 ) identity++;
				}
			}
			weight[k*nseq+l] = (overlap - identity)*1.0/maxlen;
			if  ( verbose & VERB_DEBUG )
				cout << "k=" << k << " l=" << l << " weight=" << weight[k*nseq+l] << endl;
		}
	}
	
	return weight;
}

/**
@brief 	Correlated mutation analysis of an alignment.
@param 	seqs			the set of sequences.
@param 	seqflag	 		boolean array to select sequence positions.
@param 	refseqid		reference sequence to report on.
@param 	cutoff			cutoff for reporting correlated mutations.
@param 	&simat			similarity matrix.
@return Matrix 			the analysis result matrix.

Reference: Gobel, Sander & Schneider (1994) Proteins 18, 309-317.
Mutation (residue variation) correlation is defined as:
					1
	r(i,j) =  ------------- sum(w(k,l)*(s(i,k,l) - <s(i)>)*(s(j,k,l) - <s(j)>))
			  m^2*o(i)*o(j)
	where:
		m:         number of sequences
		o(i):      standard deviation of similarities at alignment position i
		w(k,l):    weight for sequences k and l
				   (1 - fractional identity: see function seq_aligned_identity)
		s(i,k,l):  similarity for alignment position i between sequences k and l
		<s(i)>:    average similarity at alignment position i
Individual high-scoring correlations (using the given cutoff value) are reported
as follows:
	Res1	Num1	Res2	Num2	Total	Corr
	T	9	I	17	210	 0.631
	TAIIIVVVIVVVIVIIIIIII
	IILLLLLLLLLLLLLLLLLLL
The first 4 values gives the type and alignment position of the correlating residues.
The total is the number of comparisons made: maximally m*(m-1)/2
The last number is the correlation coefficient.
The following two lines gives the corresponding residues at the two alignment positions
for all the sequences, allowing the user to see on what basis this is a high correlation.

**/
Matrix		sequence_correlated_mutation(vector<Bsequence>& seqs, vector<int> seqflag,
					string& refseqid, double cutoff, Bresidue_matrix& simat)
{
	long 				h, i, j, k, l;
	long				refseqnum(0), nseq(seqs.size());
	
	for ( auto& seq: seqs ) {
		if ( seq.identifier().find(refseqid) != string::npos )
			break;
		refseqnum++;
	}

	Matrix				mat;

	if ( refseqnum >= nseq ) {
		cerr << "Sequence " << refseqid << " not found!" << endl;
		return mat;
	}
	
	long				maxlen = sequence_maximum_length(seqs);
	
	Bsequence&			seqref = seqs[refseqnum];
	
	string&				refseq = seqref.sequence();
	vector<char>		selseq(maxlen*nseq, 0);

	if  ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_correlated_mutation: maxlen=" << maxlen << " nseq=" << nseq << endl;
		
	int 				m1, m2, ir, ntot;
	int 				total = nseq*(nseq-1)/2;		// Only lower triangle
	double				sum, ssum, maxCC(-1e30);
	vector<double>		avg(maxlen, 0);
	vector<double>		std(maxlen, 0);
	vector<double>		num(maxlen, 0);
	vector<double>		score(maxlen*nseq*nseq, 0);
	
	string				code1 = simat.code();
	Matrix				sim_mat = simat.matrix();
    
	int 				igap = code1.find('-');
	
	if  ( verbose & VERB_DEBUG )
		cout << "DEBUG seq_correlated_mutation: total=" << total << endl;
	
	vector<double>		weight = sequence_aligned_weight(seqs);
	
	for ( i=ir=0; i<maxlen; i++ ) {
		if ( seqflag[i] ) {
			sum = ssum = 0;
			m1 = code1.find(seqref.sequence()[i]);
			for ( k=1; k<nseq; ++k ) {
				m1 = code1.find(seqs[k].sequence()[i]);
				for ( l=0; l<k; ++l ) {
					m2 = code1.find(seqs[l].sequence()[i]);
					h = (ir*nseq + k)*nseq + l;
						// Exclude all gaps
					if ( m1 != igap && m2 != igap ) {
						score[h] = sim_mat[m2][m1];
						sum += score[h];
						ssum += score[h]*score[h];
						num[ir]++;
					} else score[h] = -999;
//					if ( m1 == m2 ) score[h] = -999;
				}
			}
			avg[ir] = std[ir] = 0;
			if ( total ) {
				avg[ir] = sum/total;
				std[ir] = sqrt(ssum/total - avg[ir]*avg[ir]);
			}
			refseq[ir] = refseq[i];
			for ( j=0; j<nseq; ++j )
				selseq[ir*nseq+j] = seqs[j].sequence()[i];
			if  ( verbose & VERB_DEBUG )
				cout << "i=" << ir << " avg=" << avg[ir] << " std=" << std[ir] << " num=" << num[ir] << endl;
			ir++;
		}
	}
	maxlen = ir;

	if ( verbose & VERB_PROCESS ) {
		cout << "Correlated mutation analysis:" << endl;
		cout << "Matrix size:                    " << maxlen << " x " << maxlen << endl;
		cout << "Reference sequence:             " << seqref.identifier() << " (" << refseqnum + 1 << ")" << endl << endl;
		cout << endl << "Res1\tNum1\tRes2\tNum2\tTotal\tCorr" << endl;
	}

	mat = Matrix(maxlen, maxlen);
	
	ir = 0;
//	sigma2 = 0;
	for ( i=0; i<maxlen; i++ ) {
		for ( j=0; j<=i; j++ ) {
			ntot = 0;
			if ( num[i] > 2 && num[j] > 2 ) {
				if ( std[i] && std[j] ) {
					for ( k=1; k<nseq; k++ ) {
						for ( l=0; l<k; l++ ) {
							m1 = (i*nseq + k)*nseq + l;
							m2 = (j*nseq + k)*nseq + l;
							if ( score[m1] > -100 && score[m2] > -100 ) {
								mat[i][j] += weight[k*nseq+l]*
									(score[m1] - avg[i])*(score[m2] - avg[j]);
								ntot++;
							}
						}
					}

					mat[i][j] /= total*std[i]*std[j];
				}
			}
			if ( i == j ) mat[i][j] = 1;
			else if ( maxCC < mat[i][j] ) maxCC = mat[i][j];
//			sigma2 += mat[i][j] * mat[i][j]; 		// Sum of squared correlations
			if ( mat[i][j] > cutoff && i != j ) {
				cout << refseq[j] << tab << j+1 << tab << 
					refseq[i] << tab << i+1 << tab << ntot << tab << mat[i][j] << endl;
				if ( verbose & VERB_PROCESS ) {
					for ( k=0; k<nseq; k++ )
						cout << selseq[j*nseq+k];
					cout << endl;
					for ( k=0; k<nseq; k++ )
						cout << selseq[i*nseq+k];
					cout << endl;
				}
				ir++;
			}
			mat[j][i] = mat[i][j];
		}
	}
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Correlations reported:           " << ir << endl;
		cout << "Maximum off-diagonal coefficient: " << maxCC << endl << endl;
	}
	
	return mat;
}

