/**
@file	ps_views.h
@brief	Header file for postscript tools dealing with views
@author 	Bernard Heymann
@date	Created: 20011127
@date	Modified: 20230512
**/

#include <fstream>
#include "Bstring.h"
#include "View2.h"
#include "Euler.h"
#include "symmetry.h"

//Function prototypes
int 		ps_views(string filename, string symmetry_string, vector<View2<double>>& views, int flags);
int 		ps_views(ofstream* fps, string symmetry_string, vector<View2<double>>& views, int flags);
int 		ps_views(string filename, vector<View2<double>>& views);
int 		ps_views(ofstream* fps, vector<View2<double>>& views);
int			ps_sets_of_views(string filename, string title, vector<View2<double>>& views, vector<double>& fom);
int			ps_phi_theta_plot(ofstream* fps, int left, int bottom, int width, int height, int ncol);


