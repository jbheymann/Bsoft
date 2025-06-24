/**
@file	Bsequence.h
@brief	Header for handling molecular sequences..
@author	Bernard Heymann
@date	Created: 20250509
@date	Modified: 20250516
**/

#include "utilities.h"

#ifndef _Bsequence_

#define	NNEIGHBOR	64

class Bsequence {
private:
	string			id;		// Identifier
	string			dsc;	// Description
	string			tp;		// Type: DNA, RNA or Protein
	string			s;		// Sequence of characters
	double			fom;	// Sequence figure-of-merit
	long			sel;	// Selection flag
public:
	Bsequence() { }
	Bsequence(string as) : id(as), fom(1), sel(1) { }
	void			identifier(string as) { id = as; }
	string&			identifier() { return id; }
	void			description(string as) { dsc = as; }
	string&			description() { return dsc; }
	void			type(string as) { tp = as; }
	string&			type() { return tp; }
	void			sequence(string as) { s = as; }
	void			add_sequence(string as) { s += as; }
	string&			sequence() { return s; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	void			FOM(double d) { fom = d; }
	long			FOM() { return fom; }
	long			length() { return s.length(); }
	char&			operator[](long i) { return s[i]; }
};

#define _Bsequence_
#endif
