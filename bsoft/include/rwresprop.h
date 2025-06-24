/**
@file	rwresprop.h
@brief	Header file for reading residue property files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20250601
**/

#include "Bstring.h"
#include "Matrix.h"

/* Constants */
#define MAXRES	25 		// Maximum number of different residues

#ifndef _Brestype_
/************************************************************************
@Object: struct Bresidue_type
@Description:
	Residue type property structure.
@Features:
	A residue type may be identified by either a one- or three-letter code.
*************************************************************************/
class Bresidue_type {
private:
	char		c1;			// Single letter representation
	char		c3[3];		// Three-letter representation
	float		ms;			// Residue mass
	float		vl;			// Residue volume
	float		ex;			// Residue extension
	float		exs;		// Residue extension standard deviation
	float		ch; 		// Residue charge
	float		hp;			// Residue hydrophobicity
	float		comp[5];	// Composition: HCNOS
public:
	Bresidue_type() { }
	Bresidue_type(char c, double m, double v, double e, double es, double chrg, double h ) {
		c1 = c; ms = m; vl = v; ex = e; exs = es; ch = chrg; hp = h;
	}
	void		code1(char c) { c1 = c; }
	void		code1(string s) { c1 = s[0]; }
	char		code1() { return c1; }
	void		code3(string s) { c3[0] = s[0]; c3[1] = s[1]; c3[2] = s[2]; }
	string		code3() { string s("UNK"); s[0]=c3[0]; s[1]=c3[1]; s[2]=c3[2]; return s; }
	void		mass(double m) { ms = m; }
	double		mass() { return ms; }
	void		volume(double v) { vl = v; }
	double		volume() { return vl; }
	void		extension(double e) { ex = e; }
	double		extension() { return ex; }
	void		extension_stdev(double es) { exs = es; }
	double		extension_stdev() { return exs; }
	void		charge(double c) { ch = c; }
	double		charge() { return ch; }
	void		hydrophobicity(double h) { hp = h; }
	double		hydrophobicity() { return hp; }
	void		composition(vector<double> v) { for ( int i=0; i<5; ++i ) comp[i] = v[i]; }
	vector<double>	composition() {
		vector<double>	v(5);
		for ( int i=0; i<5; ++i ) v[i] = comp[i];
		return v;
	}
} ;

/************************************************************************
@Object: struct Bresidue_matrix
@Description:
	Residue relationship matrix.
@Features:
	A pairwise residue relationship matrix to encode a property such 
	as similarity.
*************************************************************************/
/* Residue relationship matrix (similarity) */
class Bresidue_matrix {
private:
    string		c;			// Symbol list for matrix (n characters)
	Matrix		m; 			// nxn Matrix (similarity)
public:
	Bresidue_matrix() { }
	Bresidue_matrix(long n) {
		c.resize(n,0);
		m = Matrix(n,n);
	}
	string&		code() { return c; }
	Matrix&		matrix() { return m; }
} ;

#define _Brestype_
#endif

/* Function prototypes */
map<char,Bresidue_type>	get_residue_properties_code1(string& filename);
map<string,Bresidue_type>	get_residue_properties_code3(string& filename);
vector<Bresidue_type>	get_residue_properties(string& filename);
int 			write_residue_properties(Bstring& filename, Bresidue_type* rt);
Bresidue_matrix	get_residue_matrix(string& filename);


