/**
@file	rwmodel_param.cpp
@brief	Library routines to read and write model dynamics parameters in STAR format
@author 	Bernard Heymann
@date	Created: 20100305
@date	Modified: 20230717
**/

#include "rwmodel_param.h"
#include "rwmodel_star.h"
#include "star.h"
#include "json.h"
#include "model_tags.h"
#include "mol_tags.h"
#include "file_util.h"
#include "string_util.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


/**
@brief 	Generates material from a model file.
@param 	*model			linked list of models.
@return Bmodparam*		new model parameter structure.
**/
Bmaterial	material_from_model(Bmodel* model)
{
	Bmaterial				material(model->identifier());
	map<string,Bcomptype>&	types = material.composition();
	material.density(RHO, DA_A3);
	
	Bmodel*				mp;
	Bcomponent*			comp;
	string				el;

    for ( mp = model; mp; mp = mp->next )
		for( comp = mp->comp; comp; comp = comp->next ) {
			el = comp->element();
			if ( types.find(el) != types.end() )
					types[el].component_count_increment();
				else
					types[el].component_count(1);
		}
	
	return material;
}


/**
@brief 	Generates a model parameter file from a model.
@param 	*model			linked list of models.
@return Bmodparam*		new model parameter structure.

	All the types in a list of models are collated into a string list of
	types for a new model parameter structure.
	The distances between component types in the matrix is set to the
	smallest such distance between components of the relevant types.
	The distance constant is set to 1.
	The distance potential type is set to 2 (soft sphere).

**/
Bmodparam	model_param_generate(Bmodel* model)
{
	Bmodparam					md;
	model_param_generate(md, model);
	return md;
}

int			model_param_generate(Bmodparam& md, Bmodel* model)
{
	int				i, j, k, l1, l2, num(0), n;
	double			d, a;
	string			id;
	Bmodel*			mp;
	Bcomptype*		c;
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	Bcomponent*		comp3;
	Vector3<float>	v1, v2;
	
	if ( verbose )
		cout << "Generating a component distance matrix from " << model->identifier() << endl << endl;
	
	map<string,Bcomptype>&		ct = md.comptype;
	ct.clear();
	
	vector<vector<Blinktype>>&	lt = md.linktype;
	lt.clear();
	
	vector<vector<vector<Bangletype>>>&	at = md.angletype;
	at.clear();
	
	vector<string>				ids;
	
	for ( num=0, mp = model; mp; mp = mp->next ) {
		for ( c = mp->type; c; c = c->next ) {
			id = c->identifier();
			if ( ct.find(id) == ct.end() ) {
				c->index(num);
				ct[id] = *c;
				ids.push_back(id);
				num++;
			}
		}
	}
	
	if ( verbose )
		cout << "Number of component types:      " << num << endl << endl;

	lt.resize(num);
	for ( i=0; i<num; ++i ) {
		lt[i].resize(num);
		for ( j=0; j<num; ++j ) {
			Blinktype&	lt1 = lt[i][j];
			lt1.identifier(ids[i],0);
			lt1.identifier(ids[j],1);
			lt1.length(1e6);
			lt1.distance(1e6);
		}
	}
	
	if ( verbose )
		cout << "Link arrays set up:             " << lt.size() << endl;

	at.resize(num);
	for ( i=0; i<num; ++i ) {
		at[i].resize(num);
		for ( j=0; j<num; ++j ) {
			at[i][j].resize(num);
			for ( k=0; k<num; ++k ) {
				Bangletype&	at1 = at[i][j][k];
				at1.identifier(ids[i],0);
				at1.identifier(ids[j],1);
				at1.identifier(ids[k],2);
				at1.angle(0);
			}
		}
	}

	if ( verbose )
		cout << "Angle arrays set up:            " << at.size() << endl;

	// Set up link lengths and distances to the shortest for each link or pair
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_param_generate: Setting up link parameters" << endl;
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->comp ) {
//		cout << "model " << m->identifier() << endl;
		for ( comp1 = mp->comp; comp1->next; comp1 = comp1->next ) {
			if ( !comp1->type() ) {
				cerr << "Error: Type not defined! model = " << mp->identifier() << ":" << comp1->identifier() << endl;
				bexit(-1);
			}
			i = comp1->type()->index();
			for ( comp2 = comp1->next; comp2; comp2 = comp2->next ) {
				if ( !comp2->type() ) {
					cerr << "Error: Type not defined! model = " << mp->identifier() << ":" << comp2->identifier() << endl;
					bexit(-1);
				}
				j = comp2->type()->index();
				d = comp1->location().distance(comp2->location());
//				cout << i << tab << j << tab << d << endl;
				Blinktype&	lt1 = lt[i][j];
				if ( comp1->find_link_exists(comp2) ) {
					cout << i << tab << j << tab << d << endl;
					if ( lt1.length() > d ) lt1.length(d);
					n++;
				} else {
					if ( lt1.distance() > d ) lt1.distance(d);
				}
			}
		}
//		cout << "model " << m->identifier() << " done" << endl;
	}

	if ( verbose )
		cout << "Link reference distances:       " << n << endl;

	// Set up angles to the smallest for each triplet of component types
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG model_param_generate: Setting up angle parameters" << endl;
	for ( n=0, mp = model; mp; mp = mp->next ) if ( mp->comp ) {
//		cout << "model " << m->identifier() << endl;
		for ( comp1 = mp->comp; comp1->next; comp1 = comp1->next ) {
			for ( l1=1; l1<comp1->link.size(); l1++ ) {
				comp2 = comp1->link[l1];
				v1 = comp2->location() - comp1->location();
				for ( l2=0; l2<l1; l2++ ) {
					comp3 = comp1->link[l2];
					v2 = comp3->location() - comp1->location();
					i = comp1->type()->index();
					j = comp2->type()->index();
					k = comp3->type()->index();
					Bangletype&	at1 = at[i][j][k];
					a = v1.angle(v2);
					if ( a > M_PI/4 && ( at1.angle() < 0.001 || at1.angle() > a ) ) {
						at1.angle(a);
						Bangletype&	at2 = at[i][k][j];
						at2.angle(a);
						n++;
					}
				}
			}
		}
	}

	if ( verbose )
		cout << "Angle references:               " << n << endl;

//	cout << "model_param_generate done" << endl;
	
	return 0;
}

/**
@brief 	Sets the component type indices from the model parameters.
@param 	*model		linked list of models.
@param 	&md			model parameter structure.
@return int			0.

	The position of a type in the distance matrix type list determines
	the type index within the distance matrix itself.
	The indices of the component types within each model is then set
	to the corresponding index for the distance matrix.
	If the type is not represented in the distance matrix, its index
	is set to -1.

**/
int			model_param_set_type_indices(Bmodel* model, Bmodparam& md)
{
	Bcomptype*		ct;
	
	if ( verbose )
		cout << "Setting type indices for a distance parameter matrix" << endl << endl;
	
	for ( ; model; model = model->next ) {
		for ( ct = model->type; ct; ct = ct->next ) {
			ct->index(-1);
			for ( auto cts: md.comptype )
				if ( cts.first == ct->identifier() ) ct->index(cts.second.index());
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_param_set_type_indices: " << ct->identifier() << " " << ct->index() << endl;
		}
	}
	
	return 0;
}

/**
@brief 	Updates model reference parameters from dynamics parameters.
@param 	*model		linked list of models.
@param 	&md			model parameter structure.
@return int			0.
**/
int			model_update_reference_parameters(Bmodel* model, Bmodparam& md)
{
	long			i, j, k, n(md.comptype.size());
	Bcomptype*		ct;
	Bcomponent*		comp;
	Blink*			link;
	Bangle*			angle;
	
	if ( verbose )
		cout << "Updating reference parameters" << endl << endl;
	
	for ( ; model; model = model->next ) {
		for ( ct = model->type; ct; ct = ct->next ) {
			ct->index(-1);
			for ( auto cts: md.comptype )
				if ( cts.first == ct->identifier() ) {
					Bcomptype&		ctr = cts.second;
					ct->index(ctr.index());
					if ( ctr.mass() ) ct->mass(ctr.mass());
					if ( ctr.charge() ) ct->charge(ctr.charge());
					if ( ctr.radius() ) ct->radius(ctr.radius());
				}
			if ( ct->index() < 0 )
				cerr << "Error: Component type " << ct->identifier() << " not found!" << endl;
			if ( verbose & VERB_DEBUG )
				cout << "DEBUG model_update_reference_parameters: " << ct->identifier() << " " << ct->index() << endl;
		}
		if ( model->link ) {
			if ( verbose )
				cerr << "Warning: Deleting original link list!" << endl;
			model->clear_links();
		}
		for ( comp = model->comp; comp; comp = comp->next ) {
			for ( auto& l: comp->link ) {
				link = model->add_link(comp, l);
				link->select(md.distancetype);
				i = comp->type()->index();
				if ( i<0 || i>=n )
					cerr << "Component type index " << i << " not found!" << endl;
				j = l->type()->index();
				if ( j<0 || j>=n )
					cerr << "Component type index " << j << " not found!" << endl;
				if ( i>=0 && j>=0 && i<n && j<n )
					link->length(md.linktype[i][j].length());
			}
		}
		if ( model->angle ) {
			if ( verbose )
				cerr << "Warning: Deleting original angle list!" << endl;
			model->clear_angles();
		}
		for ( comp = model->comp; comp; comp = comp->next ) {
			i = comp->type()->index();
			if ( i<0 || i>=n )
				cerr << "Component type index " << i << " not found!" << endl;
			for ( auto& l1: comp->link ) {
				j = l1->type()->index();
				if ( j<0 || j>=n )
					cerr << "Component type index " << j << " not found!" << endl;
				for ( auto& l2: comp->link ) {
					if ( l1->identifier() < l2->identifier() ) {
						k = l2->type()->index();
						if ( k<0 || k>=n )
							cerr << "Component type index " << k << " not found!" << endl;
						if ( i>=0 && j>=0 && k>=0 && i<n && j<n && k<n ) {
							angle = model->add_angle(comp, l1, l2);
							angle->angle(md.angletype[i][j][k].angle());
						}
					}
				}
			}
		}
	}
	
	return 0;
}

Bstar		read_parameter_file(string& filename)
{
	string			atfile;
	string			propfile;
	if ( filename.length() ) atfile = filename;
	else atfile = "atom_prop.star";
	
	if ( access(atfile.c_str(), R_OK) == 0 )
		propfile = atfile;
	else
		propfile = parameter_file_path(atfile);

	filename = propfile;

	if ( filename.length() < 1 ) {
		cerr << "Error: No parameter file specified!" << endl;
		bexit(-1);
	}
	
 	Bstar					star;
	
 	if ( star.read(filename) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 )
		cerr << "No data blocks found in the STAR file!" << endl;

	return star;
}

map<string,Bcomptype> 	read_atom_properties(Bstar& star)
{
	map<string,Bcomptype>	ct;

	if ( star.blocks().size() < 0 ) return ct;
	
	BstarBlock			block = star.block(0);

	if ( block.loops().size() < 0 ) return ct;
	
	BstarLoop			loop = block.loop(0);

	// Check if it contains atom parameters
	if ( loop.find(ATOM_TYPE_SYMBOL) < 0 ) {
		cerr <<  "Error: Star to atom properties conversion unsuccessful!" << endl;
		return ct;
	}
	
	string				symbol;
	long				j, k;
	vector<double>		coef(11,0);
	
	// Convert the atomic data
	for ( auto il: block.loops() ) {
		if ( ( j = il.find(ATOM_TYPE_SYMBOL) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				symbol = ir[j];
				Bcomptype	c(symbol);
				if ( ( k = il.find(ATOM_TYPE_NUMBER) ) >= 0 ) c.index(to_integer(ir[k]));
				if ( ( k = il.find(ATOM_TYPE_MASS) ) >= 0 ) c.mass(to_real(ir[k]));
				if ( ( k = il.find(ATOM_TYPE_OXIDATION) ) >= 0 ) c.charge(to_real(ir[k]));
				if ( ( k = il.find(ATOM_TYPE_RADIUS_BOND) ) >= 0 ) c.link_radius(to_real(ir[k]));
				if ( ( k = il.find(ATOM_TYPE_RADIUS_VDW) ) >= 0 ) c.radius(to_real(ir[k]));
				if ( ( k = il.find(ATOM_TYPE_SCAT_A1) ) >= 0 ) coef[0] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_A2) ) >= 0 ) coef[1] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_A3) ) >= 0 ) coef[2] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_A4) ) >= 0 ) coef[3] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_A5) ) >= 0 ) coef[4] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_B1) ) >= 0 ) coef[5] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_B2) ) >= 0 ) coef[6] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_B3) ) >= 0 ) coef[7] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_B4) ) >= 0 ) coef[8] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_B5) ) >= 0 ) coef[9] = to_real(ir[k]);
				if ( ( k = il.find(ATOM_TYPE_SCAT_C) ) >= 0 ) coef[10] = to_real(ir[k]);
				c.coefficients(coef);
				ct[symbol] = c;
//				cout << symbol << tab << c.index() << endl;
//				c.show();
			}
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of atom parameters read: " << ct.size() << endl;
	
	return ct;
}

map<string,Bcomptype> 	read_atom_properties(string& filename)
{
	Bstar				star = read_parameter_file(filename);

	return read_atom_properties(star);
}

map<string,Bmaterial> 	read_material_star(string& filename)
{
	Bstar				star = read_parameter_file(filename);
	
	map<string,Bmaterial>	material;

	if ( star.blocks().size() < 0 ) return material;
	
	BstarBlock			block = star.block(0);

	// Check if it contains material parameters
	if ( !block.exists(MATERIAL_NAME) ) {
		cerr << "Error: Star to material properties conversion unsuccessful! (" << filename << ")" << endl;
		return material;
	}
	
	if ( block.loops().size() < 0 ) return material;
	
	BstarLoop			loop = block.loop(0);

	// Check if it contains atom parameters
	if ( loop.find(ATOM_TYPE_SYMBOL) < 0 ) {
		cerr << "Error: Star to material properties conversion unsuccessful! (" << filename << ")" << endl;
		return material;
	}

	string 				atompropfile;		// Atom properties file
	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);

	if ( verbose & VERB_FULL )
		cout << "Converting material parameters" << endl;
	
	string				symbol;
	long				j, k, n(0);
	
	// Convert the material objects
	for ( auto ib: star.blocks() ) {
//		ib.show_tags();
		Bmaterial		m;
		string			name(ib.at(MATERIAL_NAME));
		m.identifier(name);
//		m.density(to_real(ib.at(MATERIAL_DENSITY)), 0);	// Units = g/cm3
		m.density(ib.real(MATERIAL_DENSITY), G_CM3);	// Units = g/cm3
		m.select(1);
		for ( auto il: ib.loops() ) {
//			il.show_tags();
			if ( ( j = il.find(ATOM_TYPE_SYMBOL) ) >= 0 ) {
				for ( auto ir: il.data() ) {
					symbol = ir[j];
					m[symbol].identifier(symbol);
					if ( ( k = il.find(ATOM_TYPE_COUNT) ) >= 0 ) m[symbol].component_count(to_integer(ir[k]));
//					cout << m[symbol].identifier() << tab << m[symbol].component_count() << endl;
				}
			}
		}
		m.update_parameters(atompar);
		material[name] = m;
		n++;
	}
	
	if ( verbose & VERB_FULL )
		cout << "Materials converted: " << n << endl << endl;

	if ( verbose & VERB_FULL )
		for ( auto m: material )
			cout << m.second.identifier() << tab << m.second.density(DA_A3) << endl;

	return material;
}

map<string,Bmaterial> 	read_material_json(string filename)
{
	JSvalue			js = JSparser(filename).parse();
	
	map<string,Bmaterial>	material;

	if ( js.size() < 1 ) return material;
	
//	JSvalue			jm = js[0];

	// Check if it contains material parameters
	if ( !js[0].exists(MATERIAL_NAME) ) {
		cerr << "Error: JSON to material properties conversion unsuccessful! (" << filename << ")" << endl;
		return material;
	}
	
	// Check if it contains atom parameters
	if ( !js[0].exists(ATOM_TYPE) ) {
		cerr << "Error: JSON to material properties conversion unsuccessful! (" << filename << ")" << endl;
		return material;
	}

	string 				atompropfile;		// Atom properties file
	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);

	if ( verbose )
		cout << "Converting material parameters" << endl;
	
	string				symbol;
	long				n(0);
	
	// Convert the material objects
	for ( auto jm: js.array() ) {
		Bmaterial		m;
		string			name(jm[MATERIAL_NAME].value());
		m.identifier(name);
		m.density(jm[MATERIAL_DENSITY].real(), G_CM3);	// Units = g/cm3
		m.select(1);
		for ( auto ja: jm[ATOM_TYPE].array() ) {
			if ( ja.exists(ATOM_TYPE_SYMBOL) ) {
				symbol = ja[ATOM_TYPE_SYMBOL].value();
				m[symbol].identifier(symbol);
			}
			if ( ja.exists(ATOM_TYPE_COUNT) )
				m[symbol].component_count(ja[ATOM_TYPE_COUNT].real());
		}
		m.update_parameters(atompar);
		material[name] = m;
		n++;
	}
	
	if ( verbose )
		cout << "Materials converted: " << n << endl;
	
	return material;
}

map<string,Bmaterial> 	read_material_properties(string& filename)
{
	string			ext = extension(filename);
	
	if ( verbose )
		cout << "Reading " << filename << endl << endl;
	
	if ( ext == "star" ) {
		return read_material_star(filename);
	} else if ( ext == "json" ) {
		return read_material_json(filename);
	}
	
	cerr << "Error: Material properties file extension " << ext << " not supported!" << endl;
	
	return map<string,Bmaterial>();
}

/**
@brief 	Reads material properties from multiple files.
@param 	file_list		list of file names.
@param	paramfile		atom properties file name.
@param	flags			1=addhydrogens
@return FileType			enumerated file type.

	The file extension is the main determinant of the file type.
	Flags:
		1	add hydrogens to model files.

**/
map<string,Bmaterial> 	read_material_properties(vector<string> file_list, string paramfile, int flags)
{
	FileType				ft;
	string					ext;
	map<string,Bmaterial>	material;

	for ( auto filename: file_list ) {
		ft = file_type(filename);
		ext = extension(filename);
		if ( verbose )
			cout << "Reading " << filename << endl;
		if ( ft == Material ) {
			map<string,Bmaterial>	m = read_material_properties(filename);
			material.insert(m.begin(), m.end());
		} else if ( ft == Model || ft == Molecule ) {
			Bmodel*	model = read_model(filename, paramfile);
			Bmaterial	m = material_from_model(model);
			model_kill(model);
			if ( flags & 1 ) material_protein_add_hydrogens(m);
			material[filename] = m;
		} else {
			cerr << "Error: Material properties file extension " << ext << " not supported!" << endl;
		}
	}
	
	map<string,Bcomptype>	atompar = read_atom_properties(paramfile);

	for ( auto& m: material ) {
		Bmaterial&		m1 = m.second;
		if ( m1.density(DA_A3) < 0.01 ) m1.density(RHO, DA_A3);
		m1.unit(DA_A3);					// Convert the density to Da/A3
		m1.update_parameters(atompar);
		if ( verbose & VERB_FULL )
			m1.show();
	}

	return material;
}

int			write_material_star(string filename, map<string,Bmaterial> material)
{
	string			id("Material_parameters");
 	Bstar			star;
	
	star.line_length(200);
	
	for ( auto it: material ) {
		Bmaterial&		m = it.second;
		BstarBlock&		block = star.add_block(m.identifier());
		block[MATERIAL_NAME] = m.identifier();
		block[MATERIAL_DENSITY] = to_string(m.gram_per_cm3());
		BstarLoop&		loop = block.add_loop();
		loop.tags()[ATOM_TYPE_SYMBOL] = 0;
		loop.tags()[ATOM_TYPE_COUNT] = 1;
		for ( auto ct: m.composition() ) {
			vector<string>&	vs = loop.add_row(2);
			vs[0] = ct.second.identifier();
			vs[1] = to_string(ct.second.component_count());
		}
		
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_material_star: " << filename << endl;

	return star.write(filename);
}

int			write_material_json(string filename, map<string,Bmaterial> material)
{
	JSvalue			js(JSarray);
	
	for ( auto it: material ) {
		Bmaterial&		m = it.second;
		JSvalue			jm(JSobject);
		jm[MATERIAL_NAME] = m.identifier();
		jm[MATERIAL_DENSITY] = m.gram_per_cm3();
		JSvalue			jt(JSarray);
		for ( auto ct: m.composition() ) {
			JSvalue		ja(JSobject);
			ja[ATOM_TYPE_SYMBOL] = ct.second.identifier();
			ja[ATOM_TYPE_COUNT] = ct.second.component_count();
			jt.push_back(ja);
		}
		jm[ATOM_TYPE] = jt;
		js.push_back(jm);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_material_json: " << filename << endl;

	js.write(filename);
	
	return 0;
}

int			write_material_properties(string filename, map<string,Bmaterial> material)
{
	string			ext = extension(filename);
	
	if ( verbose )
		cout << "Writing " << filename << endl << endl;
	
	if ( ext == "star" ) {
		return write_material_star(filename, material);
	} else if ( ext == "json" ) {
		return write_material_json(filename, material);
	}
	
	cerr << "Error: Material properties file extension " << ext << " not supported!" << endl;
	
	return -1;
}

/**
@brief 	Reads a model parameter file.
@param 	filename		file name.
@return Bmodparam		new model parameter structure.

	The only format supported is a STAR format.

**/
Bmodparam 	read_dynamics_parameters(string filename)
{
	Bmodparam					md;
	update_dynamics_parameters(md, filename);
	return md;
}

/**
@brief 	Updates model parameters from a parameter file.
@param 	&md				model parameters.
@param 	filename		file name.
@return int				0.

	The only format supported is a STAR format.

**/
int			update_dynamics_parameters(Bmodparam& md, string filename)
{
	// Atom parameter file
	string			atfile;
	string			propfile;
	if ( filename.length() ) atfile = filename;
	else atfile = "atom_prop.star";
	
	if ( access(atfile.c_str(), R_OK) == 0 )
		propfile = atfile;
	else
		propfile = parameter_file_path(atfile);

	filename = propfile;

	map<string,Bcomptype>&		ct = md.comptype;
	ct.clear();
	
	vector<vector<Blinktype>>&	lt = md.linktype;
	lt.clear();
	
	vector<vector<vector<Bangletype>>>&	at = md.angletype;
	at.clear();

	if ( filename.length() < 1 ) {
		cerr << "Error: No parameter file specified!" << endl;
		bexit(-1);
	}
	
 	Bstar					star;
	
 	if ( star.read(filename) < 0 ) {
		error_show(filename.c_str(), __FILE__, __LINE__);
		return -1;
	}
	
	if ( verbose )
		cout << "Reading parameter file " << filename << endl << endl;
	
	BstarBlock			block = star.block(0);
	BstarLoop			loop = block.loop(0);

	// Check if it contains atom or component parameters
//	if ( loop.find(ATOM_TYPE_SYMBOL) < 0 && loop.find(COMPTYPE_ID) < 0 ) {
//		cerr <<  "Error: Star to component properties conversion unsuccessful!" << endl;
//		return -1;
//	}

	string				atompropfile;
	long				ctmx(0);
	if ( loop.find(ATOM_TYPE_SYMBOL) < 0 && loop.find(COMPTYPE_ID) < 0 )
		ct = read_atom_properties(atompropfile);
		
	for ( auto& c: ct ) if ( ctmx < c.second.index() ) ctmx = c.second.index();

	string				symbol;
	long				n(0), i, j, k;
	vector<double>		coef(11,0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_dynamics_parameters: Reading " << filename << endl;
	
	// Convert the atomic data
	for ( auto il: block.loops() ) {
//		il.show_tags();
		if ( ( j = il.find(COMPTYPE_ID) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				symbol = ir[j];
				Bcomptype	c(symbol);
				if ( ctmx < n ) ctmx = n;
				c.index(n++);
				if ( ( k = il.find(COMPTYPE_FILENAME) ) >= 0 ) c.file_name(ir[k]);
				if ( ( k = il.find(COMPTYPE_NUMBER) ) >= 0 ) c.image_number(to_integer(ir[k]));
				if ( ( k = il.find(COMPTYPE_MASS) ) >= 0 ) c.mass(to_real(ir[k]));
				if ( ( k = il.find(COMPTYPE_CHARGE) ) >= 0 ) c.charge(to_real(ir[k]));
				if ( ( k = il.find(COMPTYPE_FOM) ) >= 0 ) c.FOM(to_real(ir[k]));
				if ( ( k = il.find(COMPTYPE_SELECT) ) >= 0 ) c.select(to_integer(ir[k]));
				ct[symbol] = c;
			}
			if ( verbose & VERB_FULL )
				cout << "Number of component types = " << ct.size() << endl;
		}

//		if ( ( j = il.find(ATOM_TYPE_SYMBOL) ) >= 0 ) {
//			ct = read_atom_properties(star);
//		}
		
		lt.resize(ctmx+1);
		for ( auto &ltr: lt ) ltr.resize(ctmx+1);

		if ( ( j = il.find(LINKTYPE_ID1) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				Blinktype	b;
				if ( ( k = il.find(LINKTYPE_ID1) ) >= 0 ) b.identifier(ir[k],0);
				if ( ( k = il.find(LINKTYPE_ID2) ) >= 0 ) b.identifier(ir[k],1);
				if ( ( k = il.find(LINKTYPE_LENGTH) ) >= 0 ) b.length(to_real(ir[k]));
				if ( ( k = il.find(LINKTYPE_DISTANCE) ) >= 0 ) b.distance(to_real(ir[k]));
				if ( ( k = il.find(LINKTYPE_KLENGTH) ) >= 0 ) b.Klength(to_real(ir[k]));
				if ( ( k = il.find(LINKTYPE_KDISTANCE) ) >= 0 ) b.Kdistance(to_real(ir[k]));
				if ( ( k = il.find(LINKTYPE_FOM) ) >= 0 ) b.FOM(to_real(ir[k]));
				if ( ( k = il.find(LINKTYPE_SELECT) ) >= 0 ) b.select(to_integer(ir[k]));
//				cout << "link " << b.identifier(0) << tab << b.identifier(1) << endl;
//				cout << "indices " << ct[b.identifier(0)].index() << tab << ct[b.identifier(1)].index() << endl;
				i = ct[b.identifier(0)].index();
				j = ct[b.identifier(1)].index();
//				cout << lt.size() << tab << lt[0].size() << endl;
				lt[i][j] = lt[j][i] = b;
			}
			if ( verbose & VERB_FULL )
				cout << "Number of link types = " << lt.size()*lt.size() << endl;
		}

		if ( ( j = il.find(BOND_TYPE_SYMBOL1) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				Blinktype	b;
				if ( ( k = il.find(BOND_TYPE_SYMBOL1) ) >= 0 ) b.identifier(ir[k],0);
				if ( ( k = il.find(BOND_TYPE_SYMBOL2) ) >= 0 ) b.identifier(ir[k],1);
				if ( ( k = il.find(BOND_TYPE_LENGTH) ) >= 0 ) b.length(to_real(ir[k]));
				if ( ( k = il.find(BOND_TYPE_VDWDIST) ) >= 0 ) b.distance(to_real(ir[k]));
				i = ct[b.identifier(0)].index();
				j = ct[b.identifier(1)].index();
				lt[i][j] = lt[j][i] = b;
			}
		}

		at.resize(ctmx+1);
		for ( auto &atr: at ) {
			atr.resize(ctmx+1);
			for ( auto &atr2: atr ) atr2.resize(ctmx+1);
		}

		if ( ( j = il.find(ANGLETYPE_ID1) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				Bangletype	a;
				if ( ( k = il.find(ANGLETYPE_ID1) ) >= 0 ) a.identifier(ir[k],0);
				if ( ( k = il.find(ANGLETYPE_ID2) ) >= 0 ) a.identifier(ir[k],1);
				if ( ( k = il.find(ANGLETYPE_ID3) ) >= 0 ) a.identifier(ir[k],2);
				if ( ( k = il.find(ANGLETYPE_ANGLE) ) >= 0 ) a.angle(to_real(ir[k])*M_PI/180.0);
				if ( ( k = il.find(ANGLETYPE_KANGLE) ) >= 0 ) a.Kangle(to_real(ir[k]));
				if ( ( k = il.find(ANGLETYPE_FOM) ) >= 0 ) a.FOM(to_real(ir[k]));
				if ( ( k = il.find(ANGLETYPE_SELECT) ) >= 0 ) a.select(to_integer(ir[k]));
				i = ct[a.identifier(0)].index();
				j = ct[a.identifier(1)].index();
				k = ct[a.identifier(2)].index();
				at[i][j][k] = at[i][k][j] = a;
			}
		}

		if ( ( j = il.find(ANGLE_TYPE_SYMBOL1) ) >= 0 ) {
			for ( auto ir: il.data() ) {
				Bangletype	a;
				if ( ( k = il.find(ANGLE_TYPE_SYMBOL1) ) >= 0 ) a.identifier(ir[k],0);
				if ( ( k = il.find(ANGLE_TYPE_SYMBOL2) ) >= 0 ) a.identifier(ir[k],1);
				if ( ( k = il.find(ANGLE_TYPE_SYMBOL3) ) >= 0 ) a.identifier(ir[k],2);
				if ( ( k = il.find(ANGLE_TYPE_ANGLE) ) >= 0 ) a.angle(to_real(ir[k])*M_PI/180.0);
				i = ct[a.identifier(0)].index();
				j = ct[a.identifier(1)].index();
				k = ct[a.identifier(2)].index();
				at[i][j][k] = at[i][k][j] = a;
			}
		}
	}
	
	if ( verbose ) {
		cout << "Model parameters:               " << filename << endl;
		cout << "Component types:                " << md.comptype.size() << endl;
		cout << "Link types:                     " << md.linktype.size() << endl;
		cout << "Angle types:                    " << md.angletype.size() << endl;
	}

	return 0;
}


/**
@brief 	Writes a model parameter file.
@param 	filename	file name.
@param 	&md			model parameter structure.
@return int			0.

The only format supported is a STAR format.

**/
int		 	write_dynamics_parameters(string filename, Bmodparam& md)
{
	string			id("Model_parameters");
 	Bstar			star;

	md.comment += command_line_time2();
	star.comment(md.comment);

//	star.comment(model->comment().str());
	star.line_length(200);
	
	BstarBlock&		block = star.add_block(id);

	if ( md.comptype.size() ) {
		BstarLoop&			loop = block.add_loop();
		loop.tags() = comptype_tags();
		for ( auto ct: md.comptype ) {
			vector<string>&	vs = loop.add_row();
			vs[0] = ct.second.identifier();
			vs[1] = ct.second.file_name();
			vs[2] = to_string(ct.second.image_number());
			vs[3] = to_string(ct.second.mass());
			vs[4] = to_string(ct.second.charge());
			vs[5] = to_string(ct.second.FOM());
			vs[6] = to_string(ct.second.select());
		}
	}
	
	if ( md.linktype.size() ) {
		BstarLoop&			loop = block.add_loop();
		loop.tags() = linktype_tags();
		for ( auto r: md.linktype ) {
			for ( auto lt: r ) if ( lt.length() < 20 || lt.distance() < 20 ) {
				vector<string>&	vs = loop.add_row();
				vs[0] = lt.identifier(0);
				vs[1] = lt.identifier(1);
				vs[2] = to_string(lt.length());
				vs[3] = to_string(lt.distance());
				vs[4] = to_string(lt.Klength());
				vs[5] = to_string(lt.Kdistance());
				vs[6] = to_string(lt.FOM());
				vs[7] = to_string(lt.select());
			}
		}
	}

	if ( md.angletype.size() ) {
		BstarLoop&			loop = block.add_loop();
		loop.tags() = angletype_tags();
		for ( auto r1: md.angletype ) {
			for ( auto r2: r1 ) {
				for ( auto at: r2 ) if ( at.angle() > 0 ) {
					vector<string>&	vs = loop.add_row();
					vs[0] = at.identifier(0);
					vs[1] = at.identifier(1);
					vs[2] = at.identifier(2);
					vs[3] = to_string(at.angle()*180.0/M_PI);
					vs[4] = to_string(at.Kangle());
					vs[5] = to_string(at.FOM());
					vs[6] = to_string(at.select());
				}
			}
		}
	}

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_star: " << filename << endl;

	int		err = star.write(filename);
	
	if ( err < 0 ) return err;
	
	return 0;
}


