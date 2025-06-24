/**
@file	ps_sequence.h 
@brief	Header file for postscript output for sequence analysis functions.
@author 	Bernard Heymann 
@date	Created: 20010515
@date	Modified: 20250602
**/
 
#include "Bmodel.h"
#include "Complex.h"
#include "utilities.h"
/*
class Bdomain {
private:
	long		st, en;
	RGBA<float> col;
public:
	Bdomain() { }
	Bdomain(long s, long e, RGBA<float> c) : st(s), en(e), col(c) { }
	void			start(long i) { st = i; }
	long			start() { return st; }
	void			end(long i) { en = i; }
	long			end() { return en; }
	void			color(RGBA<float> c) { col = c; }
	RGBA<float>&	color() { return col; }
} ;
*/

// Function prototypes
int 		ps_seq_info(Bstring& filename, Bstring& title, int nres, int length,
				double* info, double* nseq, Complex<float>* per, double* freq, char* pattern);
int 		ps_seq_info(Bstring& filename, Bstring& title, int nres,
				vector<double>& info, vector<double>& nseq,
				vector<Complex<float>>& per, vector<double>& freq, string& pattern);
int 		ps_seq_hydrophob(Bstring& filename, Bstring& title, int length,
				double* Hphob, int* HPseg, double* nseq, Complex<float>* per);
int 		ps_seq_hydrophob(Bstring& filename, Bstring& title,
				vector<double>& Hphob, vector<int>& HPseg,
				vector<double>& nseq, vector<Complex<float>>& per);
int			ps_plot_domains(string& filename, string& name, long bar_width, 
				long bar_height, vector<Bgroup>& doms, bool numbers);

		
