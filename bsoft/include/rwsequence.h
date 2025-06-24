/**
@file	rwsequence.h
@brief	Header file for reading molecular sequence files
@author Bernard Heymann
@date	Created: 19980822
@date	Modified: 20250510
**/

#include "Bsequence.h"

/* Constants */
#define MAXSEQLEN	1000000	// Maximum sequence length

/* Function prototypes */
vector<Bsequence>	read_sequence(string& filename);
vector<Bsequence>	read_sequence(vector<string>& file_list);
long		write_sequence(string& filename, vector<Bsequence> seqs);
long		write_sequence(char *filename, vector<Bsequence> seqs);


/*
Bmolgroup*  molgroup_init();
Bmolecule*	molecule_add(Bmolecule** mol, char* name);
Bmolecule*	molecule_add(Bmolecule** mol, Bstring& name);
Bresidue*	residue_add(Bresidue** res, const char* type);
Bresidue*	residue_add(Bresidue** res, Bstring& type);
Batom*		atom_add(Batom** atom, const char* type);
Batom*		atom_add(Batom** atom, Bstring& type);
Batom*		atom_copy(Batom* atom);
long		residue_count(Bmolgroup* molgroup);
long		atom_count(Bmolgroup* molgroup);
int			atom_clean_type(Batom* atom, const char* type);
Bbond*		bond_add(Bbond** bond, Batom* atom1, Batom* atom2, double l, double k);
//Bangle*		angle_add(Bangle** angle, Batom* atom1, Batom* atom2, Batom* atom3, double a, double k);
int 		molgroup_list_kill(Bmolgroup* molgroup);
int 		molgroup_kill(Bmolgroup* molgroup);
int 		molecule_kill(Bmolecule* mol);
int 		residue_kill(Bresidue* res);
int			bond_kill(Bbond* bond);
//int			angle_kill(Bangle* angle);
Bmolgroup*	molgroup_list_copy(Bmolgroup* molgroup);
Bmolgroup*	molgroup_copy(Bmolgroup* molgroup);
Bmolecule* 	molecule_copy(Bmolecule* mol);
Bmolecule* 	mol_copy_and_add_to_molgroup(Bmolgroup* molgroup, Bmolecule* mol);
Bbond*		molgroup_bond_list_copy(Bmolgroup* molgroup, Bmolgroup* molgroupcopy);
int			molgroup_from_molgroup_list(Bmolgroup* molgroup);
Bmolgroup*  read_molecule(const char *filename, const char *select, const char* paramfile);
Bmolgroup*  read_molecule(Bstring& filename, Bstring& atom_select, Bstring& paramfile);
Bmolgroup*	read_molecule(Bstring* file_list, int set_pbc, Vector3<double> box, 
				Bstring atom_select, Bstring paramfile);
int 		write_molecule(char *filename, Bmolgroup* molgroup);
int 		write_molecule(Bstring& filename, Bmolgroup *molgroup);
int			molgroup_list_write(Bstring& filename, Bmolgroup* molgroup);
long 		molgroup_count_molecules(Bmolgroup* molgroup);
long 		molgroup_count_residues(Bmolgroup* molgroup);
long 		mol_count_residues(Bmolecule *mol);
long 		molgroup_count_atoms(Bmolgroup* molgroup);
long 		mol_count_atoms(Bmolecule *mol);
int			molgroup_consolidate_gaps(Bmolgroup* molgroup);
long 		molgroup_stats(Bmolgroup* molgroup, int show);
long 		molgroup_stats(Bmolgroup* molgroup);
long 		mol_stats(Bmolecule* mol, int show);
long 		mol_stats(Bmolecule* mol);
int 		molecule_update_comment(Bmolgroup* molgroup, int n, char** strings);
int 		molecule_get_masses(Bmolgroup* molgroup, Bstring& paramfile);
Bbond*		molgroup_bond_list_generate(Bmolgroup* molgroup, double maxlength, int wrap);
Bbond*		mol_bond_list_generate(Bmolgroup* molgroup, double bondlength, int wrap);
int			molecules_to_molgroups(Bmolgroup* molgroup);
*/
