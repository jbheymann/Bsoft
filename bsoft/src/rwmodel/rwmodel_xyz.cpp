/**
@file	rwmodel_xyz.cpp
@brief	Library routines to read and write Kirkland's xyz atomic model parameters
@author 	Bernard Heymann
@date	Created: 20220426
@date	Modified: 20221115
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "model_mol.h"
#include "string_util.h"
#include "utilities.h"
#include <fstream>

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads Kirkland's xyz atomic model parameters.
@param 	*file_list		list of model parameter file names.
@param 	&atompar		parameters for atomic Z numbers.
@return Bmodel*			model parameters.
**/
Bmodel*		read_model_xyz(vector<string> file_list, map<string,Bcomptype>& atompar)
{
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Bcomptype*		ct = NULL;
	string			id("1"), type("XYZ");
	ifstream		fmod;
	string			s;
	long			Z, natom(0);
	double			occ, wobble; // Occupancy and Debye-Waller (thermal motion)
	Vector3<double>	size, loc;

//	map<string,Bcomptype> 	atompar = read_atom_properties(paramfile);
	map<long,string>		Zsymbol;
	
	for ( auto a: atompar ) Zsymbol[a.second.index()] = a.first;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG read_model_xyz: Number of atoms: " << atompar.size() << endl;
		for ( auto z: Zsymbol )
			cout << z.first << tab << z.second << endl;
	}
	
	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return  NULL;
		if ( model ) mp = model->add(base(filename));
		else mp = model = new Bmodel(base(filename));
		mp->model_type(type);
		mp->select(1);
		getline(fmod, s);	// Description
		if ( verbose & VERB_PROCESS )
			cout << "Description:  " << s << endl;
		fmod >> size[0] >> size[1] >> size[2];
		while ( fmod.good() ) {
			fmod >> Z >> loc[0] >> loc[1] >> loc[2] >> occ >> wobble;
			if ( fmod.good() && Z > 0 ) {
				if ( Zsymbol.find(Z) == Zsymbol.end() ) {
					cerr << "Atom with number " << Z << " not found!" << endl;
					bexit(-1);
				}
				natom++;
				if ( comp ) comp = comp->add(natom);
				else comp = model->comp = new Bcomponent(natom);
				comp->location(loc);
				comp->density(occ);
				comp->FOM(wobble);
				comp->select(1);
				ct = model->add_type(Zsymbol[Z]);
				comp->type(ct);
				ct->index(Z);
				comp->description(Zsymbol[Z]);
			}
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes Kirkland's xyz atomic model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_xyz(string& filename, Bmodel* model, int splt)
{
	int					n, Z(1);
	Bmodel*				mp = NULL;
	Bcomponent*			comp;
	string				onename;
	double				occ(1), wobble(0);
	Vector3<double>		loc, size, off;
	string				paramfile;

	map<string,Bcomptype> 	atompar = read_atom_properties(paramfile);

	model->calculate_bounds();
	size = model->maximum() - model->minimum();
	if ( size[1] > size[0] ) size[0] = size[1];
	else size[1] = size[0];
//	size[0] = size[1] = 1.5*size[0];
//	off = (size - model->maximum() - model->minimum())*0.5;
	
	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
//			onename = filename.pre_rev('.') + string(n+1, format) + filename.post_rev('.');
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << command_line().c_str() << endl;
		fmod << setw(16) << size[0] << setw(16) << size[1] << setw(16) << size[2] << endl;
 		for ( comp = mp->comp; comp; comp = comp->next ) {
 			Z = atompar[comp->element()].index();
			loc = comp->location() + off;
//			occ = comp->density();
//			wobble = comp->FOM();
			fmod << setw(5) << Z << setw(14) << loc[0] << setw(14)
				<< loc[1] << setw(14) << loc[2] << setw(14) << occ
				<< setw(14) << wobble << endl;
		}
		fmod << setw(5) << -1 << endl;
		fmod.close();
	}
	
	return 0;
}
