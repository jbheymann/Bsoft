/**
@file	rwmodel_gro.cpp
@brief	Library routines to read and write GROMACS coordinate files
@author 	Bernard Heymann
@date	Created: 19980822
@date	Modified: 20250318
**/

#include "rwmodel_gro.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Gromacs coordinate files.
@param 	file_list		coordinate file name list.
@param 	&atompar		parameters.
@return Bmodel*			model parameters.

	Gromacs format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	MD of 2 waters, t= 0.0
	    6
	    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
	    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
	    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180

**/
Bmodel*		read_model_gro(vector<string> file_list, map<string,Bcomptype>& atompar)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_gro: filename=" << file_list[0] << endl;
	    
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	string			id("1"), type("GRO");
	ifstream		fmod;
	string			s, el, restype, atomtype;
	long			n(0), natom(0), atomnum, resnum;
	Vector3<double>	loc, vel;
	
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
		getline(fmod, s);	// Description
		if ( verbose & VERB_PROCESS )
			cout << "Description:  " << s << endl;
		getline(fmod, s);	// Atom
		natom = to_integer(s);
		if ( verbose & VERB_PROCESS )
			cout << "Number of atoms:  " << natom << endl;
		while ( !fmod.eof() ) {
			getline(fmod, s);	// Atom
			if ( s.size() < 40 ) break;
			resnum = to_integer(s.substr(0, 5));
			restype = s.substr(5, 5);
			restype = remove_spaces(restype);
			atomtype = s.substr(10, 5);
			atomtype = remove_spaces(atomtype);
			atomnum = to_integer(s.substr(15, 5));
			loc[0] = to_real(s.substr(20, 8));
			loc[1] = to_real(s.substr(28, 8));
			loc[2] = to_real(s.substr(36, 8));
			if ( s.size() > 40 ) {
				vel[0] = to_real(s.substr(44, 8));
				vel[1] = to_real(s.substr(52, 8));
				vel[2] = to_real(s.substr(60, 8));
			}
			natom++;
			if ( comp ) comp = comp->add(atomnum);
			else comp = mp->comp = new Bcomponent(atomnum);
			comp->location(loc*10);
			comp->velocity(vel*10);
			comp->density(1);
			comp->FOM(1);
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
@brief 	Writes a Gromacs coordinate file.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int 			number of models written (<0 if writing failed).

	Gromacs format:
	01234567890123456789012345678901234567890123456789012345678901234567890
	MD of 2 waters, t= 0.0
	    6
	    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
	    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
	    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180

**/
int 	write_model_gro(string& filename, Bmodel* model, int splt)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_gro: filename=" << filename << endl;
	
	if ( !model ) return -1;
	    
	long				n(0), nc;
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
		nc = mp->component_count();
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << command_line().c_str() << endl;
		fmod << nc << endl;
 		for ( comp = mp->comp; comp; comp = comp->next ) {
 			vector<string>&	vd = comp->description();
 			loc = comp->location()/10;
 			vel = comp->velocity()/10;
			fmod << setw(5) << vd[4] << left << setw(5) << vd[2] << right << setw(5)
				<< vd[1] << setw(5) << comp->identifier() << fixed << setprecision(3)
				<< setw(8) << loc[0] << setw(8) << loc[1] << setw(8) << loc[2] 
				<< setw(8) << vel[0] << setw(8) << vel[1] << setw(8) << vel[2] << endl;
		}
		fmod.close();
	}
	
	return n;
}

