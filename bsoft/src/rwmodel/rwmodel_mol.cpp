/**
@file	rwmodel_mol.cpp
@brief	Library routines to read and write atomic model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20250514
**/

#include "rwmodel.h"
//#include "rwmolecule.h"
//#include "model_links.h"
//#include "model_util.h"
#include <fstream>
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads MDL molfile model parameters.
@param 	*file_list		list of model parameter file names.
@param 	&atompar		parameters for atomic Z numbers.
@return Bmodel*			model parameters.
**/
Bmodel*		read_model_mol(vector<string> file_list, map<string,Bcomptype>& atompar)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_mol: filename=" << file_list[0] << endl;
	    
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;
//	Bcomptype*		ct = NULL;
//	string			id("1"), type("WAH");
	ifstream		fmod;
	string			s, el;
	long			atomnum, natom(0), nlink(0);
	long			x, y, z, l1, l2;
	map<long,Bcomponent*>	htab;
	
	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		fmod.open(filename.c_str());
		if ( fmod.fail() ) return  NULL;
		getline(fmod, s);	// Title
		if ( model ) mp = model->add(s);
		else mp = model = new Bmodel(s);
//		mp->model_type(type);
		mp->select(1);
		getline(fmod, s);	// Viewer?
		getline(fmod, s);	// Blank?
 		getline(fmod, s);	// Numbers?
		istringstream	ss(s);
		ss >> natom >> nlink;
		cout << "Components: " << natom << " Links: " << nlink << endl;
		atomnum = 0;
		while ( !fmod.eof() ) {
			getline(fmod, s);
			if ( s.find("END") != string::npos ) break;
			istringstream	ss(s);
			if ( atomnum < natom ) {	// Atom
				ss >> x >> y >> z >> el;
				atomnum++;
				if ( mp->comp ) comp = comp->add(atomnum);
				else comp = mp->comp = new Bcomponent(atomnum);
				comp->description(el);
				comp->location(x,y,z);
				htab[atomnum] = comp;
			} else {					// Link
				ss >> l1 >> l2;
				if ( mp->link ) link = link->add(htab[l1], htab[l2]);
				else link = mp->link = new Blink(htab[l1], htab[l2]);
			}
		}
		fmod.close();
	}

	return model;
}

/**
@brief 	Writes MDL molfile model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@param 	splt		flag to split into separate models.
@return int			models written.
**/
int			write_model_mol(string& filename, Bmodel* model, int splt)
{
	int				n;
	Bmodel*			mp = NULL;
	Bcomponent*		comp;
	Blink*			link;
	string			onename;

	ofstream		fmod;

	for ( n=0, mp = model; mp; mp = mp->next, n++ ) {
		if ( model->next )
			onename = insert(filename, n+1, splt);
		else
			onename = filename;
		fmod.open(onename.c_str());
		if ( fmod.fail() ) return  -1;
		fmod << mp->identifier() << endl;
		fmod << "  Bsoft" << endl;
		fmod << mp->comment() << endl;
		fmod << mp->component_count() << mp->link_count() << endl;
		for ( comp = mp->comp; comp; comp = comp->next )
			fmod << fixed << setprecision(4) << right <<
				setw(10) << comp->location()[0] << 
				setw(10) << comp->location()[1] << 
				setw(10) << comp->location()[2] <<
				" " << comp->description()[0] << endl;
		for ( link = mp->link; link; link = link->next )
			fmod << setw(3) << link->comp[0]->identifier() <<
				setw(3) << link->comp[0]->identifier() << endl;
		fmod << "M  END" << endl;
		fmod << endl;
		fmod.close();
	}
	
	return 0;
}

/**
@brief 	Reads molecular model parameters.
@param 	*file_list	list of model parameter file names.
@param 	&paramfile	parameter file.
@return Bmodel*		model parameters.
**/
/*Bmodel*		read_model_molecule(vector<string> file_list, string& paramfile)
{
    string    		atom_select("all");
	string			id;
	Bmolgroup*		molgroup;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;
	Bbond*			bond;
	
	int				i(0);
	RGBA<float>		rgba(1,1,1,1);
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_molecule: " << file_list[0] << endl;
	
	for ( auto filename: file_list ) {
		if ( verbose & VERB_LABEL )
			cout << "Reading file:                   " << filename << endl;
		molgroup = read_molecule(filename.c_str(), atom_select.c_str(), paramfile.c_str());
		if ( mp ) mp = mp->add(++i);
		else mp = model = new Bmodel(++i);
		mp->model_type(mp->identifier());
		if ( molgroup ) {
			if ( molgroup->id.length() && molgroup->id[0] != ' ') mp->identifier(molgroup->id.str());
			mp->symmetry(molgroup->pointgroup.str());
			comp = NULL;
			link = NULL;
			for ( bond = molgroup->bond; bond; bond = bond->next ) {
				link = new Blink();
				if ( !mp->link ) mp->link = link;
				link->radius(1);
				link->select(1);
				link->color(rgba);
			}
			for ( mol = molgroup->mol; mol; mol = mol->next ) {
				for ( res = mol->res; res; res = res->next ) {
					for ( atom = res->atom; atom; atom = atom->next ) {
						if ( comp ) comp = comp->add(atom->num);
						else comp = mp->comp = new Bcomponent(atom->num);
						comp->location(atom->coord);
						comp->FOM(atom->b);
						comp->select((int) (atom->q + 0.999));
						id = atom->type;
//						id = id.no_space();
						comp->type(mp->add_type(id));
						for ( bond = molgroup->bond, link = mp->link; bond && link; bond = bond->next, link = link->next ) {
							if ( atom == bond->atom1 ) link->comp[0] = comp;
							else if ( atom == bond->atom2 ) link->comp[1] = comp;
						}
					}
				}
			}
			molgroup_kill(molgroup);
		}
	}

	models_setup_links(model);
	
	return model;
}
*/

/**
@brief 	Writes molecular model parameters.
@param 	&filename	model parameter file name.
@param 	*model		model parameters.
@return int			models written.
**/
/*int			write_model_molecule(string& filename, Bmodel* model)
{
	int				i, n;
	Bstring			restype("UNK");
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		mg = NULL;
	Bmolgroup*		mgc = NULL;
	Bmolecule*		mol = NULL;
	Bresidue*		res = NULL;
	Batom*  		atom = NULL, *atom2;
	Bbond*			bond = NULL;
	
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Blink*			link = NULL;

	for ( n=i=0, mp = model; mp; mp = mp->next, n++ ) {
		mg = molgroup_init();
		if ( !molgroup ) molgroup = mg;
		else mgc->next = mg;
		mgc = mg;
//		mg->comment = model->comment().pre('\n') + "\nREMARK Model";
		mg->id = mp->identifier();
		mg->pointgroup = mp->symmetry();
		mol = molecule_add(&mg->mol, mg->id);
		res = residue_add(&mol->res, restype);
		res->num = 1;
		atom = NULL;
		for ( comp = mp->comp; comp; comp = comp->next ) {
			atom = atom_add(&atom, comp->type()->identifier().c_str());
			if ( !res->atom ) res->atom = atom;
			res->insert[0] = ' ';
			atom->num = stoi(comp->identifier());
			atom->coord = comp->location();
			atom->sel = comp->select();
			atom->q = comp->select();
			atom->b = comp->FOM();
		}
		bond = NULL;
		for ( link = mp->link; link; link = link->next ) {
			for ( comp = mp->comp, atom = res->atom; comp && comp != link->comp[0] && atom; 
				comp = comp->next, atom = atom->next ) ;
			for ( comp = mp->comp, atom2 = res->atom; comp && comp != link->comp[1] && atom2; 
				comp = comp->next, atom2 = atom2->next ) ;
			if ( atom && atom2 ) {
				bond = (Bbond *) bond_add(&bond, atom, atom2, link->length(), 1);
				if ( !mg->bond ) mg->bond = bond;
			} else {
				cerr << "Error: Problem with link " << link->comp[0]->identifier() << " to " << link->comp[1]->identifier() << endl;
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_molecule: atoms converted: " << i << endl;

	Bstring			fn(filename);
	molgroup_list_write(fn, molgroup);

	molgroup_list_kill(molgroup);
	
	return  n;
}
*/
