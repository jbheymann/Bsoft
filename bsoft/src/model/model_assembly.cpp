/**
@file	model_assembly.cpp
@brief	Library routines used for model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250318
**/

#include "model_assembly.h"
#include "model_transform.h"
//#include "mol_transform.h"
//#include "mol_compare.h"
//#include "mol_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
long		molgroup_write_into_grid(Bmolgroup* molgroup, Vector3<int> size, Vector3<double> min, double sampling, int* grid)
{
	long			i, x, y, z;
	Bmolecule*		mol;
	Bresidue*		res;
	Batom*			atom;

	long	vol = (long) size.volume();
	char*			tgrid = new char[vol];
	for ( i=0; i<vol; i++ ) tgrid[i] = 0;
	
	for ( mol=molgroup->mol; mol; mol=mol->next ) {
		for ( res=mol->res; res; res=res->next ) {
			for ( atom=res->atom; atom; atom=atom->next ) {
				x = (int) ((atom->coord[0] - min[0])/sampling);
				y = (int) ((atom->coord[1] - min[1])/sampling);
				z = (int) ((atom->coord[2] - min[2])/sampling);
				i = (z*size[1]+y)*size[0]+x;
				if ( i < 0 || i > vol ) {
					cerr << "Error in molgroup_write_into_grid: i=" << i << " (vol=" << vol << ")" << endl;
					return -1;
				}
				tgrid[i] = 1;
			}
		}
	}
	
	for ( i=0; i<vol; i++ ) grid[i] += tgrid[i];

	delete[] tgrid;
	
	return 0;
}
*/

/**
@brief 	Concatenates selected molecules into one group.
@param 	*model		model parameters.
@param 	paramfile	atomic parameter file.
@param 	separate	flag to generate separate molecule groups.
@return Bmolgroup*	list of molecule groups.

	Only the first model in the linked list is processed.

**/
/*Bmolgroup*	model_assemble(Bmodel* model, string paramfile, int separate)
{
	Bcomponent*		comp = NULL;
//	Bcomptype*		comptype = NULL;

	Bmolgroup*		mglist = NULL;
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		molgroup1 = NULL;
	Bmolecule*		mol = NULL;
    string    		atom_select("all");
	Quaternion		q;
	Transform		t;

	long			i, nsel, nover;
	double			sampling = 5;
	Vector3<double>	min(model->comp->location()), max(model->comp->location());
	Vector3<int>	size;
	
	if ( verbose )
		cout << "Assembling components" << endl;
	
	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		min = min.min(comp->location());
		max = max.max(comp->location());
		nsel++;
	}
	
	string 		fn(model->type->file_name());
	molgroup1 = read_molecule(fn.c_str(), atom_select.c_str(), paramfile.c_str());
	min -= molgroup1->box;
	max += molgroup1->box;
	molgroup_kill(molgroup1);
	for ( i=0; i<3; i++ ) size[i] = (int) ((max[i] - min[i])/sampling);
	
	if ( nsel < 1 ) {
		cerr << "Error: No components are selected!" << endl;
		return NULL;
	}
	
	long	vol = (long) size.volume();
	int*			grid = new int[vol];
	for ( i=0; i<vol; i++ ) grid[i] = 0;
	
	for ( nsel=0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		nsel++;
//		comptype = model_get_type(model, comp->type);
		fn = comp->type()->file_name();
		molgroup1 = read_molecule(fn.c_str(), atom_select.c_str(), paramfile.c_str());
		molgroup1->id = comp->identifier();
//		q = quaternion_from_view(comp->view);
		q = comp->view().quaternion();
//		t = transform_from_quaternion(q);
		t = Transform(q);
		t.origin = molgroup_center_of_mass(molgroup1);
		t.trans = comp->location() - t.origin;
		molgroup_coor_rotate(molgroup1, t);
		molgroup_stats(molgroup1);
		if ( molgroup_write_into_grid(molgroup1, size, min, sampling, grid) < 0 ) {
			error_show("Error in model_assemble", __FILE__, __LINE__);
			return NULL;
		}
		if ( separate ) {
			if ( !mglist ) mglist = molgroup = molgroup1;
			else {
				molgroup->next = molgroup1;
				molgroup = molgroup1;
			}
		} else {
			if ( molgroup ) {
				if ( mol ) {
					for ( ; mol->next; mol = mol->next ) ;
					mol->next = molgroup1->mol;
				} else {
					mol = molgroup->mol = molgroup1->mol;
				}
				molgroup1->mol = NULL;
				molgroup_kill(molgroup1);
			} else {
				mglist = molgroup = molgroup1;
				mol = molgroup->mol;
			}
		}
	}

	for ( i=nover=0; i<vol; i++ ) if ( grid[i] > 1 ) nover++;
	
	delete[] grid;
		
	if ( verbose ) {
		cout << "Components assembled:           " << nsel << endl;
		cout << "Molecule group overlap:         " << 
			sampling*sampling*sampling*nover << " A3 (" << nover*100.0/size.volume() << " %)" << endl << endl;
	}
	
	return mglist;
}
*/

Bmodel*		model_assemble(Bmodel* model, string paramfile)
{
	Bmodel*			numod = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
	Quaternion		q;
	Transform		t;
	long			nsel(0);
	string			fn;
	
	if ( verbose )
		cout << "Assembling components" << endl;
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		nsel++;
		fn = comp->type()->file_name();
		mp = read_model(fn, paramfile);
		if ( numod ) numod->add(mp);
		else numod = mp;
		q = comp->view().quaternion();
		t = Transform(q);
		t.origin = mp->center_of_coordinates();
		t.trans = comp->location() - t.origin;
		model_rotate(mp, t);
		mp->calculate_bounds();
	}

	if ( verbose ) {
		cout << "Components assembled:           " << nsel << endl;
//		cout << "Molecule group overlap:         " << 
//			sampling*sampling*sampling*nover << " A3 (" << nover*100.0/size.volume() << " %)" << endl << endl;
	}
	
	return numod;
}

/**
@brief 	Calculates the centers-of-mass of molecule group components and generates a new model.
@param 	*molgroup	list of molecule groups.
@return Bmodel*		new model.

	Each molecule is assumed to be a component.

**/
/*Bmodel*		model_generate_com(Bmolgroup* molgroup)
{
	string			id, path;
	string			comptype("VER");
	Bmolgroup*		mg;
	Bmolecule*		mol;
	
	int				i, j, n=0;
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Generating a centers-of-mass model" << endl << endl;
	
	for ( i=1, mg = molgroup; mg; mg = mg->next, i++ ) {
//		if ( mg->id.length() ) mp->identifier(mg->id.str());
//		else mp->identifier() = to_string(i);
		if ( mg->id.length() ) id = mg->id.str();
		else id = to_string(i);
		if ( model ) mp = mp->add(id);
		else mp = model = new Bmodel(id);
		comp = NULL;
		for ( j=1, mol = molgroup->mol; mol; mol = mol->next, j++, n++ ) {
			cout << "Adding molecule " << j << " as component" << endl;
//			comp = component_add(&comp, j);
//			if ( !mp->comp ) mp->comp = comp;
			if ( comp ) comp = comp->add(j);
			else mp->comp = comp = new Bcomponent(j);
			comp->location(mol_center_of_mass(mol));
			if ( mol->id.length() ) id = mol->id.no_space().str();
			else id = comptype;
//			comp->type = model_add_type_by_id_and_filename(mp, id, molgroup->filename, 0);
			comp->type(mp->add_type(id, molgroup->filename.c_str(), 0));
		}
	}
	
	cout << "Models generated:               " << --i << endl;
	cout << "Components generated:           " << n << endl << endl;

	model_check(model, path);
	
	return model;
}
*/

/**
@brief 	Generates an assembly from a set of models.
@param	file_list	list of file names corresponding to the models.
@param 	paramfile	parameter file.
@return Bmodel*		new model.

	The center of coordinates is calculated from each model.
	A component is generated from each model.

**/
Bmodel*		model_generate_assembly(vector<string> file_list, string paramfile)
{
	string			id("Assembly"), path;
	long			i(0);
	Bmodel*			model = new Bmodel(id);
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;

	if ( verbose & VERB_PROCESS )
		cout << "Generating an assembly model" << endl << endl;
	
	for ( auto fn: file_list ) {
		i++;
		mp = read_model(fn, paramfile);
		if ( mp->identifier().length() ) id = mp->identifier();
		else id = to_string(i);
		if ( verbose )
			cout << "Adding model " << id << " as component" << endl;
		if ( comp ) comp = comp->add(id);
		else model->comp = comp = new Bcomponent(id);
		comp->location(mp->center_of_coordinates());
		comp->type(mp->add_type(id, fn, 0));
	}
	
	if ( verbose )
		cout << "Models assembled:               " << i << endl;

	model_check(model, path);
	
	return model;
}

/**
@brief 	Finds the molecule views with respect to a reference.
@param 	*model		model parameters.
@param 	&reffile	reference molecule file name.
@param 	&paramfile	atomic parameter file.
@return long			number of molecules selected.

	The positioning of each molecule is based on the center of mass of the reference.

**/
/*long		model_find_views(Bmodel* model, Bstring& reffile, Bstring& paramfile)
{
	Bcomponent*		comp = NULL;
	Transform		t;
    Bstring    		atom_select("all");
	
	if ( !model->type ) {
		cerr << "Error: No component types found!" << endl;
		return -1;
	}
	
	long			nsel(0);
	Bstring			fn;
	Bmolgroup*		molgroup = NULL;
	Bmolgroup*		ref_molgroup = read_molecule(reffile, atom_select, paramfile);
	
	if ( verbose )
		cout << "Molecule\tDist\t\t\tAxis\t\t\tAngle" << endl;
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() > 0 ) {
		nsel++;
		fn = comp->type()->file_name();
		molgroup = read_molecule(fn, atom_select, paramfile);
		t = molgroup_find_transformation(molgroup, ref_molgroup);
		comp->view(View2<float>(t.angle, t.axis));
		molgroup_kill(molgroup);
		if ( verbose )
			cout << comp->type()->file_name()  << tab << comp->location() << tab << comp->view() << endl;
	}
	
	molgroup_kill(ref_molgroup);
	
	return nsel;
}
*/
