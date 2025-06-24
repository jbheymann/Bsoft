/**
@file	rwmodel_cif.cpp
@brief	Library routines to read and write molecule files in CIF format
@author 	Bernard Heymann
@date	Created: 19980822 
@date	Modified: 20250619
**/

#include "star.h"
#include "rwmolecule.h"
#include "rwmodel_cif.h"
#include "mol_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Reads a molecule group from a CIF format file.
@param 	file_list		a list of file names.
@return long				number of molecules read (<0 if reading failed).
**/
Bmodel*		read_model_cif(vector<string> file_list)
{
 	Bstar			star;
	
	if ( verbose )
		cout << "Reading file:                   " << file_list[0] << endl;

	for ( auto f: file_list ) star.read(f);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_model_cif: blocks: " << star.blocks().size() << endl;
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return 0;
	}

	long			i, j, nmol(0), nres(0), natom(0), rs, re;
	string			entity, dsc, pdb, typestr, sid, gid, chain, el, atomtype;
	string			restype, restype2, resnum, resinsert;
	map<string,int>	sheet_strands;
	map<string,int>	sheet_order;

	Bmodel*			model_list = NULL;
	Bmodel*			model = NULL;
	Bcomponent*		comp = NULL;
//	Bcomponent*		comp2 = NULL;
//	Blink*			link = NULL;
	Bgroup*			group = NULL;

	for ( auto ib: star.blocks() ) {
//		model->identifier(ib.tag());
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_model_cif: block tag: " << ib.tag() << endl;
/*		molname = ib.tag();
		if ( molname.length() < 1 ) molname = "A";
		mol = molecule_add(&molgroup->mol, molname);
		mol->seq = ib.at(MOLECULE_SEQUENCE).c_str();
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_model_cif: sequence " << mol->seq << endl;
		mol->nres = ib.real(MOLECULE_LENGTH);
		if ( mol->seq.length() > 0 )
			if ( mol->nres < 1 ) mol->nres = mol->seq.length();*/

		// Make sure a new model is created
		if ( ib.exists(MOL_ID) ) {
			chain = ib.at(MOL_ID);
			if ( model_list ) model = model_list->add(chain);
			else model_list = model = new Bmodel(chain);
			if ( ib.exists(MOL_ENTITY) )
				model->model_type(ib.at(MOL_ENTITY));
			if ( ib.exists(ENTITY_ID) )
				entity = ib.at(ENTITY_ID);
			if ( ib.exists(ENTITY_DESCRIPTION) )
				dsc = ib.at(ENTITY_DESCRIPTION);
			model->description(dsc);
			if ( verbose & VERB_DEBUG ) {
				cout << "DEBUG read_model_cif: Molecule " << model->identifier() << " Type: " << model->model_type() << endl;
				cout << "DEBUG read_model_cif: Description -" << model->description() << "-" << endl;
			}
		} else for ( auto il: ib.loops() ) {
			if ( il.find(MOL_ID) >= 0 ) {
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(MOL_ID) ) >= 0 )
						chain = ir[i];
					if ( model_list ) model = model_list->add(chain);
					else model_list = model = new Bmodel(chain);
					if ( ( i = il.find(MOL_ENTITY) ) >= 0 )
						model->model_type(ir[i]);
					if ( verbose & VERB_DEBUG )
						cout << "DEBUG read_model_cif: Molecule " << model->identifier() << " Type: " << model->model_type() << endl;
				}
			}
		}
		if ( !model )
			cerr << "Error: No model found!" << endl;
		// Fill in the model information
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_model_cif: models read so far: " << model_list->count() << endl;
		for ( auto il: ib.loops() ) {
			if ( il.find(ENTITY_ID) >= 0 ) {
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(ENTITY_ID) ) >= 0 )
						entity = ir[i];
					if ( ( i = il.find(ENTITY_DESCRIPTION) ) >= 0 )
						dsc = ir[i];
					for ( model = model_list; model; model = model->next )
						if ( model->model_type() == entity ) {
							model->description(dsc);
							if ( verbose & VERB_PROCESS )
								cout << model->identifier() << tab << entity << tab << model->description() << endl;
						}
				}
/*			} else if ( il.find(CHEMICAL_ATOM_ID) >= 0 ) {
				chain = ib.at(CHEMICAL_ATOM_ID);
				if ( model ) model = model->add(chain);
				else model_list = model = new Bmodel(chain);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG read_model_cif: Chemical " << model->identifier() << endl;
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(CHEMICAL_ATOM_NUMBER) ) >= 0 ) {
						if ( comp ) comp = comp->add(ir[i]);
						else comp = model->comp = new Bcomponent(ir[i]);
						comp->density(1);
						comp->select(1);
						natom++;
					}
					if ( ( i = il.find(CHEMICAL_ATOM_SYMBOL) ) >= 0 )
						el = ir[i];
					if ( ( i = il.find(CHEMICAL_ATOM_ID) ) >= 0 )
						atomtype = ir[i];
					if ( ( i = il.find(CHEMICAL_ATOM_RES) ) >= 0 )
						restype = ir[i];
//					typestr = el + " " + atomtype + " " + restype;
					comp->type(model->add_type(atomtype));
					comp->description(el);
					comp->add_description(atomtype);
					comp->add_description(restype);
					if ( ( i = il.find(CHEMICAL_ATOM_X) ) >= 0 )
						comp->location()[0] = to_real(ir[i]);
					if ( ( i = il.find(CHEMICAL_ATOM_Y) ) >= 0 )
						comp->location()[1] = to_real(ir[i]);
					if ( ( i = il.find(CHEMICAL_ATOM_Z) ) >= 0 )
						comp->location()[2] = to_real(ir[i]);
					if ( ( i = il.find(CHEMICAL_ATOM_CHARGE) ) >= 0 )
						comp->charge(to_real(ir[i]));
				}*/
			}
			if ( il.find(HELIX_ID) >= 0 ) {
//				il.show_tags();
				j = 0;
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(HELIX_ID) ) >= 0 )
						gid = to_string(++j) + " " + ir[i];
					if ( ( i = il.find(HELIX_CHAIN) ) >= 0 )
						chain = ir[i];
					if ( ( i = il.find(HELIX_RESIDUE1) ) >= 0 )
						restype = ir[i];
					if ( ( i = il.find(HELIX_RESIDUE2) ) >= 0 )
						restype2 = ir[i];
					if ( ( i = il.find(HELIX_RESNUM1) ) >= 0 )
						rs = to_integer(ir[i]);
					if ( ( i = il.find(HELIX_RESNUM2) ) >= 0 )
						re = to_integer(ir[i]);
					model = model_list->find(chain);
					group = model->add_group(gid);
					group->group_type("HELIX ");
					group->type("1");								// Type of helix
					// Chain residue1 residue2
					group->description(chain + " " + restype + " " + restype2);
//					seq = resnum + "-" + resnum2;					// Residue numbers
//					group->sequence(seq);
					group->start(rs);							// Residue numbers
					group->end(re);
				}
			}
			if ( il.find(SHEET_ORDER_ID) >= 0 ) {
//				il.show_tags();
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(SHEET_ORDER_ID) ) >= 0 )
						sid = ir[i];
					if ( ( i = il.find(SHEET_ORDER_STRAND_ID) ) >= 0 )
						gid = ir[i] + " " + sid;
					if ( ( i = il.find(SHEET_ORDER_SENSE) ) >= 0 ) {
						if ( ir[i].find("anti-parallel") != string::npos ) sheet_order[gid] = -1;
						else sheet_order[gid] = 1;
					}
//					cout << gid << tab << ir[i] << tab << sheet_order[gid] << endl;
				}
			}
			if ( il.find(SHEET_ID) >= 0 ) {
//				il.show_tags();
				for ( auto ir: il.data() ) {
					if ( ( i = il.find(SHEET_ID) ) >= 0 )
						sid = ir[i];
					if ( ( i = il.find(SHEET_STRAND_ID) ) >= 0 )
						gid = ir[i] + " " + sid;
					if ( ( i = il.find(SHEET_CHAIN) ) >= 0 )
						chain = ir[i];
					if ( ( i = il.find(SHEET_RESIDUE1) ) >= 0 )
						restype = ir[i];
					if ( ( i = il.find(SHEET_RESIDUE2) ) >= 0 )
						restype2 = ir[i];
					if ( ( i = il.find(SHEET_RESNUM1) ) >= 0 )
						rs = to_integer(ir[i]);
					if ( ( i = il.find(SHEET_RESNUM2) ) >= 0 )
						re = to_integer(ir[i]);
					model = model_list->find(chain);
					group = model->add_group(gid);
					group->group_type("SHEET ");
					group->type(to_string(sheet_order[gid]));		// Direction of strand
					// Chain residue1 residue2 strands_in_sheet
					group->description(chain + " " + restype + " " + restype + " ");
//					seq = resnum + "-" + resnum2;					// Residue numbers
//					group->sequence(seq);
					group->start(rs);							// Residue numbers
					group->end(re);
					if ( sheet_strands.find(sid) != sheet_strands.end() )
						sheet_strands[sid] += 1;
					else
						sheet_strands[sid] = 1;
				}
				for ( model = model_list; model; model = model->next ) {
					for ( group = model->group; group; group = group->next ) {
						if ( group->group_type() == "SHEET " ) {
							sid = group->identifier();
							sid = sid.substr(sid.find(" ")+1);
//							cout << sid << tab << sheet_strands[sid] << endl;
							group->description() += to_string(sheet_strands[sid]);
						}
					}
				}
			}
			if ( il.find(ATOM_NUMBER) >= 0 ) {
//				il.show_tags();
				if ( ( i = il.find(ATOM_CHAIN) ) >= 0 )
					chain = il.data()[0][i];
				model = model_list->find(chain);
				if ( !model )
					cerr << "Error: No chain " << chain << " found!" << endl;
				comp = model->comp;
				for ( auto ir: il.data() ) {
//					for ( auto it: ir ) cout << it << endl;
					if ( ( i = il.find(ATOM_CHAIN) ) >= 0 ) {
						if ( chain != ir[i] ) {
							chain = ir[i];
							model = model_list->find(chain);
							if ( !model )
								cerr << "Error: No chain " << chain << " found!" << endl;
							comp = model->comp;
						}
					}
					if ( ( i = il.find(ATOM_NUMBER) ) >= 0 ) {
						if ( comp ) comp = comp->add(ir[i]);
						else comp = model->comp = new Bcomponent(ir[i]);
						comp->density(1);
						comp->select(1);
						natom++;
					}
					if ( ( i = il.find(ATOM_PDB) ) >= 0 )
						if ( ir[i] == "HETATM" ) comp->select(2);
					if ( ( i = il.find(ATOM_ENTITY) ) >= 0 ) {
						long	si = to_integer(ir[i]);
						if ( si ) comp->select(si);
					}
//					if ( ( i = il.find(ATOM_PDB) ) >= 0 )
//						cout << "-" << ir[i] << "-" << endl;
					if ( ( i = il.find(ATOM_SYMBOL) ) >= 0 )
						el = ir[i];
					if ( ( i = il.find(ATOM_TOPTYPE) ) >= 0 )
						atomtype = ir[i];
					if ( ( i = il.find(ATOM_RESNUMBER) ) >= 0 )
						resnum = ir[i];
					if ( ( i = il.find(ATOM_RESINSERT) ) >= 0 ) {
						resinsert = ir[i];
						if ( resinsert.size() < 1 ) resinsert = "?";
					}
					if ( ( i = il.find(ATOM_RESIDUE) ) >= 0 )
						restype = ir[i];
					comp->type(model->add_type(atomtype));
					comp->description(el);
					comp->add_description(atomtype);
					comp->add_description(restype);
					comp->add_description(chain);
					comp->add_description(resnum);
					comp->add_description(resinsert);
					if ( ( i = il.find(ATOM_X) ) >= 0 )
						comp->location()[0] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Y) ) >= 0 )
						comp->location()[1] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_Z) ) >= 0 )
						comp->location()[2] = to_real(ir[i]);
					if ( ( i = il.find(ATOM_OCCUPANCY) ) >= 0 )
						comp->density(to_real(ir[i]));
					if ( ( i = il.find(ATOM_BFACTOR) ) >= 0 )
						comp->FOM(to_real(ir[i]));
					if ( ( i = il.find(ATOM_CHARGE) ) >= 0 )
						comp->charge(to_real(ir[i]));
				}
			}
		}
//		nmol++;
//		molname[0]++;
//		if ( molname[0] > 'Z' ) molname[0] = 'A';
	}
		
	if ( nmol < 1 ) nmol--;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG read_model_cif: " << nmol << " molecules" << endl;
		cout << "DEBUG read_model_cif: " << nres << " residues" << endl;
		cout << "DEBUG read_model_cif: " << natom << " atoms" << endl;
	}
	
	return model_list;
}

long 		write_model_cif_one(string& filename, Bmodel* model, int splt)
{
 	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_cif_one: filename=" << filename << endl;

	Bstar			star;

	star.line_length(120);
	
	Bmodel*			mp = NULL;
	Bcomponent*		comp = NULL;
//	Bcomponent*		comp2 = NULL;
//	Blink*			link = NULL;
//	Bgroup*			group = NULL;

	map<string,string>	entities;

	BstarBlock&		block = star.add_block(model->identifier());

	BstarLoop&		loop1 = block.add_loop();
	loop1.tags()[MOL_ID] = 0;
	loop1.tags()[MOL_ENTITY] = 1;
	for ( mp = model; mp; mp = mp->next ) {
		vector<string>&	vs = loop1.add_row(2);
		vs[0] = mp->identifier();
		vs[1] = mp->model_type();
		entities[mp->model_type()] = mp->description();
//		cout << vs[0] << tab << vs[1] << endl;
	}
	
	BstarLoop&		loop2 = block.add_loop();
	loop2.tags()[ENTITY_ID] = 0;
	loop2.tags()[ENTITY_TYPE] = 1;
	loop2.tags()[ENTITY_DESCRIPTION] = 2;
	for ( auto e1: entities ) {
		vector<string>&	vs = loop2.add_row(3);
		vs[0] = e1.first;
		vs[1] = "polymer";
		vs[2] = '"' + e1.second + '"';
		if ( verbose & VERB_FULL ) 
			cout << vs[0] << tab << vs[2] << endl;
	}
	
	BstarLoop&		loop = block.add_loop();
	loop.tags()[ATOM_NUMBER] = 0;
	loop.tags()[ATOM_SYMBOL] = 1;
	loop.tags()[ATOM_TOPTYPE] = 2;
	loop.tags()[ATOM_RESIDUE] = 3;
	loop.tags()[ATOM_RESNUMBER] = 4;
	loop.tags()[ATOM_RESINSERT] = 5;
	loop.tags()[ATOM_CHAIN] = 6;
	loop.tags()[ATOM_ENTITY] = 7;
	loop.tags()[ATOM_X] = 8;
	loop.tags()[ATOM_Y] = 9;
	loop.tags()[ATOM_Z] = 10;
	loop.tags()[ATOM_OCCUPANCY] = 11;
	loop.tags()[ATOM_BFACTOR] = 12;
	loop.tags()[ATOM_CHARGE] = 13;

	for ( mp = model; mp; mp = mp->next ) {
//		cout << mp->identifier() << endl; 
		for ( comp = mp->comp; comp; comp = comp->next ) {
			vector<string>&	vd = comp->description();
			vector<string>&	vs = loop.add_row(14);
			vs[0] = comp->identifier();
			vs[1] = vd[0];
			vs[2] = vd[1];
			if ( vd.size() > 2 ) vs[3] = vd[2];
			else vs[3] = "UNK";
			if ( vd.size() > 4 ) vs[4] = vd[4];
			else vs[4] = to_string(1);
			if ( vd.size() > 5 ) vs[5] = vd[5];
			else vs[5] = "?";
			if ( vd.size() > 3 ) vs[6] = vd[3];
			else vs[6] = "A";
			vs[7] = mp->model_type();
			vs[8] = to_string(comp->location()[0]);
			vs[9] = to_string(comp->location()[1]);
			vs[10] = to_string(comp->location()[2]);
			vs[11] = to_string(comp->density());
			vs[12] = to_string(comp->FOM());
			vs[13] = to_string(comp->charge());
		}
	}
			
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_cif_one: " << filename << endl;

	star.write(filename, splt);
	
	return model->count();
}

long 		write_model_cif_split(string& filename, Bmodel* model, int splt)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG write_model_cif_split: filename=" << filename << endl;

	for ( Bmodel* mp = model; mp; mp = mp->next ) {
		Bmodel*		mpc = mp->copy();
		write_model_cif_one(filename, mpc, splt);
		delete mpc;
	}
		
	return model->count();
}

/**
@brief 	Writes a set of models to a CIF format file.
@param 	&filename	the file name.
@param 	*model		linked list of models.
@param 	splt		flag to split on models.
@return long 			number of models written (<0 if writing failed).
**/
long 		write_model_cif(string& filename, Bmodel* model, int splt)
{
	if ( splt ) return write_model_cif_split(filename, model, splt);
	
	return write_model_cif_one(filename, model, splt);
}

