/**
@file	rwsymop.h
@brief	Header file for reading and writing symmetry operators
@author Bernard Heymann
@date	Created: 19990509
@date	Modified: 20230623
**/

#include "symmetry.h"

// Function prototypes in symmetry.c and rwsymop.c
vector<float>	read_symat(string& filename, int spacegroup);
char*		read_symop(string& filename, int spacegroup, int& nsym);
vector<string>	read_symop(string& filename, int spacegroup);
int 		write_symat(string& filename, int spacegroup);
int 		write_pointgroup(string filename, string& symmetry_string, View2<double> ref_view);
int 		write_pointgroup(string filename, Bsymmetry& sym, View2<double> ref_view);

