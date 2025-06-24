/**
@file	rwmodel_wah.cpp
@brief	Library routines to read and write Wayne Hendrickson coordinate files
@author 	Bernard Heymann
@date	Created: 20050217
@date	Modified: 20250318
**/

#include "rwmodel_wah.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a Wayne Hendrickson coordinate file.
@param 	file_list		coordinate file name list.
@param 	&atompar		parameters.
@return Bmodel*			model parameters.

	WAH format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	   SER 1071CA     8.38771  32.21584 115.00745   0.00000   0.00000
	   ARG 1072C     11.26685  34.89560 115.61476   0.00000   0.00000

**/
Bmodel*		read_model_wah(vector<string> file_list, map<string,Bcomptype>& atompar)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_wah: filename=" << file_list[0] << endl;
	    
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	string			id("1"), type("WAH");
	ifstream		fmod;
	string			s, el, restype, atomtype;
	long			n(0), atomnum, resnum;
	Vector3<double>	loc;
	
	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return  NULL;
		if ( model ) mp = model->add(base(filename));
		else mp = model = new Bmodel(base(filename));
		n++;
		mp->model_type(type);
		mp->select(1);
		atomnum = 0;
		while ( !fmod.eof() ) {
			getline(fmod, s);	// Atom
			if ( s.size() < 60 ) break;
			restype = s.substr(3, 4);
			restype = remove_spaces(restype);
			resnum = to_integer(s.substr(7, 4));
			atomtype = s.substr(11, 4);
			atomtype = remove_spaces(atomtype);
			atomnum++;
			loc[0] = to_real(s.substr(15, 10));
			loc[1] = to_real(s.substr(25, 10));
			loc[2] = to_real(s.substr(35, 10));
			if ( comp ) comp = comp->add(atomnum);
			else comp = mp->comp = new Bcomponent(atomnum);
			comp->density(to_real(s.substr(45, 10)));
			comp->FOM(to_real(s.substr(55, 10)));
			comp->location(loc);
			comp->select(1);
			ct = mp->add_type(atomtype);
			comp->type(ct);
			el = atomtype.substr(0,1);
			comp->description(el);
			comp->add_description(atomtype);
			comp->add_description(restype);
			comp->add_description(to_string(n));
			comp->add_description(to_string(resnum));
//			component_element(comp, atompar);
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes a Wayne Hendrickson coordinate file.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int 			number of models written (<0 if writing failed).

	WAH format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	   SER 1071CA     8.38771  32.21584 115.00745   0.00000   0.00000
	   ARG 1072C     11.26685  34.89560 115.61476   0.00000   0.00000

**/
int 	write_model_wah(string& filename, Bmodel* model, int splt)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_wah: filename=" << filename << endl;
		
	if ( !model ) return -1;
	    
	long				n(0);
	Bmodel*				mp = NULL;
	Bcomponent*			comp;
	string				onename;
	Vector3<double>		loc, vel;
	string				paramfile;

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, ++n ) {
		if ( model->next )
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		for ( comp = mp->comp; comp; comp = comp->next ) {
 			vector<string>&	vd = comp->description();
 			loc = comp->location();
 			vel = comp->velocity();
			fmod << "   " << left << setw(4) << vd[2] << right << setw(4) << vd[4] 
				<< left << setw(4) << vd[1] << fixed << setprecision(5) << right
				<< setw(10) << loc[0] << setw(10) << loc[1] << setw(10) << loc[2] 
				<< setw(10) << comp->density() << setw(10) << comp->FOM() << endl;
		}
		fmod.close();
	}
	
	return n;
}

