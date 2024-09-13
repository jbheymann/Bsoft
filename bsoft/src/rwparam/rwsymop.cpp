/**
@file	rwsymop.cpp
@brief	Library routines to read and write symmetry operators
@author Bernard Heymann
@date	Created: 19991225
@date	Modified: 20230623
**/

#include "rwsymop.h"
#include "star.h"
#include "sym_tags.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Internal function prototypes
vector<string>	read_symop_star(string& filename, int spacegroup);
vector<string>	read_symop_lib(string& filename, int spacegroup);
int 		write_symop_star(string& filename, int spacegroup, vector<string>& symop, int line_len);
int 		write_pointgroup_star(string filename, Bsymmetry& sym, View2<double> ref_view);
vector<float>	sym_matrices_from_text_list(vector<string>& symop, int line_len);

/**
@brief 	Reading crystallographic symmetry operators.
@param 	&filename		file name.
@param 	spacegroup		crystal space group number.
@param 	&nsym			number of symmetry operators.
@return vector<float>		set of 12-value symmetry matrices.

	The symmetry operators are encoded as a set of matrices.

**/
vector<float>	read_symat(string& filename, int spacegroup, int& nsym)
{
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symat: Getting symmetry matrices from: " << filename << endl;
	
	vector<string>	symop = read_symop(filename, spacegroup);
	
	nsym = symop.size();
	
	// Set up the symmetry operator matrix
	vector<float>	mat = sym_matrices_from_text_list(symop, 80);
	
	return mat;
}

/**
@brief 	Reading crystallographic symmetry operators.
@param 	&symopfile		file name.
@param 	spacegroup		crystal space group number.
@param	&nsym			number of symmetry operators.
@return char*			80 character symmetry operators.

	The symmetry operators are encoded as 80 character lines.

**/
char*		read_symop(string& symopfile, int spacegroup, int& nsym)
{
	vector<string>	op = read_symop(symopfile, spacegroup);
	
	nsym = op.size();
	
	int				i, j, k, len;
	char*			symop = new char[80*nsym];
	
	for ( i=k=0; i<nsym; ++i ) {
		len = op[i].length();
		for ( j=0; j<80; ++j, ++k ) {
			if ( j<len ) symop[k] = op[i][j];
			else symop[k] = ' ';
		}
	}
	
	return symop;
}

/**
@brief 	Reading crystallographic symmetry operators.
@param 	&symopfile		file name.
@param 	spacegroup		crystal space group number.
@return vector<string>	80 character symmetry operators.

	The symmetry operators are encoded as 80 character lines.

**/
vector<string>	read_symop(string& symopfile, int spacegroup)
{
	vector<string>	symop;
	
	// No space group has more than 192 operators
	if ( spacegroup < 1 || spacegroup > 100000 ) return symop;
	
	if ( symopfile.length() < 1 ) {
		symopfile = "symop.star";
		symopfile = parameter_file_path(symopfile);
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop: Reading symmetry operator file " << symopfile << endl;

	string		ext = extension(symopfile);
	
	if ( spacegroup > 1 ) {
		if ( ext.find("star") != string::npos || ext.find("cif") != string::npos ) {
			symop = read_symop_star(symopfile, spacegroup);
		} else {
			symop = read_symop_lib(symopfile, spacegroup);
		}
		if ( !symop.size() )
			cerr << "Error: No symmetry operator file read!" << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop: Symmetry operator file read: " << symopfile << endl;
	
	return symop;
}

/**
@brief 	Writing crystallographic symmetry operators.
@param 	&filename		file name.
@param 	spacegroup		crystal space group number.
@return int				error code (<0 means failure).
**/
int 		write_symat(string& filename, int spacegroup)
{
	int				err(0);
	string			temp;
	vector<string>	symop = read_symop(temp, spacegroup);
	if ( symop.size() < 1 ) return -1;
	
	err = write_symop_star(filename, spacegroup, symop, 80);
	
	return err;
}

/**
@brief 	Writing point group symmetry operators.
@param 	&filename			file name.
@param 	&symmetry_string	symmetry string.
@param 	ref_view			reference view.
@return int					error code (<0 means failure).
**/
int 		write_pointgroup(string filename, string& symmetry_string, View2<double> ref_view)
{
	Bsymmetry	sym(symmetry_string);
	
	int			err = write_pointgroup(filename, sym, ref_view);

	return err;
}

int 		write_pointgroup(string filename, Bsymmetry& sym, View2<double> ref_view)
{
	int			err(0);

	if ( filename.find(".star") != string::npos )
		err = write_pointgroup_star(filename, sym, ref_view);
	else {
		cerr << "Error: File type for " << filename << " not supported!" << endl;
		err = -1;
	}

	return err;
}

// Find space group label and operators in a STAR format file
vector<string>	read_symop_star(string& filename, int spacegroup)
{
 	Bstar			star;
	vector<string>	symop;
	
 	if ( star.read(filename) < 0 )
		error_show(filename.c_str(), __FILE__, __LINE__);
	
	if ( star.blocks().size() < 0 ) {
		cerr << "No data blocks found in the STAR file!" << endl;
		return symop;
	}

	int				j;

	for ( auto ib: star.blocks() ) {
		if ( spacegroup == ib.integer(SYMMETRY_NUMBER) ) {
			for ( auto il: ib.loops() ) {
				if ( ( j = il.find(SYMMETRY_EQUIVXYZ) ) >= 0 ) {
//					nsym = il.data().size();
					for ( auto ir: il.data() ) {
//						string	symstr(ir[j]);
//						strcpy(symop+80*i, symstr.c_str());
						symop.push_back(ir[j]);
					}
				}
			}
			break;
		}
	}
	
	return symop;
}

// Find space group label and operators in "symop.lib" or similar file
vector<string>	read_symop_lib(string& filename, int spacegroup)
{
	int				i, j, k(0), notfound(1), number(0), nlines, nsym(0);
	char			aline[80], symall[4000];
	for ( i=0; i<4000; i++ ) symall[i] = 0;

	vector<string>	symop;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Symmetry operator library: " << filename << endl << endl;

	string			atfile;
	string			symopfile;
	if ( filename.empty() ) symopfile = filename;
	else symopfile = "symop.lib";
	
	if ( access(symopfile.c_str(), R_OK) != 0 )
		symopfile = parameter_file_path(symopfile);

	ifstream		fsym;
	fsym.open(symopfile.c_str());
	if ( fsym.fail() ) return symop;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Symmetry operator library: " << symopfile << endl << endl;
	
	while ( notfound && fsym.getline(aline, 80) ) {
		sscanf( aline, "%d %d", &number, &nlines );
		if ( number == spacegroup ) notfound = 0;
	}
	
	if ( notfound ) {
		fsym.close();
		return symop;
	}
	
	for ( i=nsym=0; i<nlines; i++ ) {
		fsym.getline(aline, 80);
		for ( j=0; j<(int)strlen(aline); j++ ) {
			if ( aline[j] != ' ' && aline[j] != '\n' ) {
				symall[k] = tolower(aline[j]);
				k++;
			}
			if ( aline[j] == '*' ) nsym++;
		}
		symall[k] = '*';
		k++;
		nsym++;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG read_symop_lib: Spacegroup=" << number << " Noperators=" << nsym << endl;

//	char*		symop = new char[nsym*80];
//	for ( i=0; i<nsym*80; i++ ) symop[i] = 0;

	j = 0;
	for ( i=0; i<nsym; i++ ) {
		k = 0;
		symop.push_back(string(80,' '));
		while( symall[j] != '*' && k<80 ) {
			symop[i][k] = symall[j];
			j++;
			k++;
		}
		j++;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG read_symop_lib: " << symop[i] << endl;
	}
	
	fsym.close();
	
	return symop;
}

int 		write_symop_star(string& filename, int spacegroup, vector<string>& symop, int line_len)
{
	int				i;
	Bstar			star;

	star.comment("# Symmetry operators\n\n");

	BstarBlock&		block = star.add_block(to_string(spacegroup));

	block[SYMMETRY_NUMBER] = to_string(spacegroup);

	BstarLoop&		loop = block.add_loop();
	loop.tags()[SYMMETRY_EQUIVID] = 0;
	loop.tags()[SYMMETRY_EQUIVXYZ] = 1;

	for ( i=0; i<symop.size(); ++i ) {
		vector<string>&	vs = loop.add_row(2);
		vs[0] = to_string(i+1);
		vs[1] = symop[i];
	}
	
	return star.write(filename);
}

int 		write_pointgroup_star(string filename, Bsymmetry& sym, View2<double> ref_view)
{
 	Bstar			star;

	star.comment("# Symmetry from bsym\n\n");

	int				i;
	Matrix3			mat = ref_view.matrix();
	Vector3<double>	v;

	BstarBlock&		block = star.add_block(sym.label());

	block[SYMMETRY_POINT_GROUP] = sym.label();
	block[SYMMETRY_PG_NUMBER] = to_string(sym.point());

	BstarLoop&		loop = block.add_loop();
	loop.tags()[SYMMETRY_AXIS_ORDER] = 0;
	loop.tags()[SYMMETRY_AXIS_X] = 1;
	loop.tags()[SYMMETRY_AXIS_Y] = 2;
	loop.tags()[SYMMETRY_AXIS_Z] = 3;

	for ( i=0; i<sym.operations(); i++ ) {
		v = mat * sym[i].axis();
		vector<string>&	vs = loop.add_row(4);
		vs[0] = to_string(sym[i].order());
		vs[1] = to_string(v[0]);
		vs[2] = to_string(v[1]);
		vs[3] = to_string(v[2]);
	}
	
	return star.write(filename);
}

/**
@brief 	Calculates symmetry matrices from a list of strings.
@param 	&symop			array of symmetry operator lines.
@param 	line_len		length of text line in the array.
@return vector<float>		a set of 12-value symmetry matrices.

	The list of strings is expected to be packed into a single character
	array with a fixed length for each string. Each string encodes a
	symmetry operation in terms of x, y and z operations in reciprocal
	space.

**/
vector<float>	sym_matrices_from_text_list(vector<string>& symop, int line_len)
{
	// Set up the symmetry operator matrix
	int 		i, j, k, l;
	int			nsym(symop.size());
	vector<float>	mat(nsym*12,0);

	char		op[200];
	for ( i=0; i<nsym; i++ ) {
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG sym_matrices_from_text_list: Symmetry operator " << i+1 << ":" << endl;
		k = 0;
		for ( j=0; j<3; j++ ) {
			l = 0;
			memset(op, 0, line_len);
			while ( k<line_len && symop[i][k] != ',' ) {
				op[l] = tolower(symop[i][k]);
				k++;
				l++;
			}
			op[l] = 0;
			if ( strstr(op,"-x") ) mat[i*12+4*j] = -1;
			else if ( strstr(op,"x") ) mat[i*12+4*j] = 1;
			if ( strstr(op,"-y") ) mat[i*12+4*j+1] = -1;
			else if ( strstr(op,"y") ) mat[i*12+4*j+1] = 1;
			if ( strstr(op,"-z") ) mat[i*12+4*j+2] = -1;
			else if ( strstr(op,"z") ) mat[i*12+4*j+2] = 1;
			if ( strstr(op,"1/2") ) mat[i*12+4*j+3] = 0.5;
			if ( strstr(op,"1/4") ) mat[i*12+4*j+3] = 0.25;
			if ( strstr(op,"3/4") ) mat[i*12+4*j+3] = 0.75;
			if ( strstr(op,"1/3") ) mat[i*12+4*j+3] = 1.0/3.0;
			if ( strstr(op,"2/3") ) mat[i*12+4*j+3] = 2.0/3.0;
			if ( strstr(op,"1/6") ) mat[i*12+4*j+3] = 1.0/6.0;
			if ( strstr(op,"5/6") ) mat[i*12+4*j+3] = 5.0/6.0;
			k++;
			if ( verbose & VERB_DEBUG )
				cout << "|" << mat[i*12+4*j] << " " << mat[i*12+4*j+1]
					<< " " << mat[i*12+4*j+2] << "|   |" << mat[i*12+4*j+3] << "|" << endl;
		}
	}
	
	return mat;
}
