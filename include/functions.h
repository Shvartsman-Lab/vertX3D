#pragma once

//DECLARATIONS
#include <stdio.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include "freeid.h"

using namespace std;

//TEMP FOR Hexagonal.cpp
int input_Nx;

//GEOMETRIC ELEMENTS
size_t Nv, Nv_pass, Ne, Nf, Nc;
double **v, **v_pass;
int **e, **f;
int **basal_edges;
int **basal_vertices, **basal_facets, **apical_facets;

//ATTRIBUTES
//vertices
double **v_F;
double **v_normal_vector;
double *v_height;
int *v_type    ;
int *v_partner ;
int *v_cell1   ;
int *v_cell2   ;
int *v_cell3   ;
int *v_cell4   ;
int *v_T1dir   ;
int *v_vertT1  ;
int *v_edgeT1  ;
double *v_clock;

//passive vertices
int *v_pass_type;
int *v_pass_cell;
int *v_pass_edge;
    
//edges
int *e_v_pass  ;
int *e_lateral1;
int *e_lateral2;
int *e_lateral3;
int *e_lateral4;
int *e_cell1   ;
int *e_cell2   ;
double *e_length;

//facets
int *f_type;
int *f_cell;
int *f_edge;
double *f_T;
    
//cells
int *c_cent_basal;
int *c_cent_apical;
double *c_V0;
double *c_kV;
double *c_alpha;
double *c_beta;

//Yolk
double VY0, c_kY;

//MISC
double *perioXYZ;
double h, Time, max_move;
double wA, wV, wY, A_tot;
FreeId v_freeId, v_pass_freeId, e_freeId, f_freeId, c_freeId;
size_t array_max, seed;
double tmax;
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
#include "allocate.h"
#include "rndom.h"
#include "torus.h"
#include "distances.h"
#include "vert_edg_fac.h"
#include "basal_network.h"
#include "basal_side.h"
#include "copy_to_apical.h"
#include "apical_side.h"
#include "lateral_side.h"
#include "output.h"
#include "initial_structure.h"
#include "force_area.h"
#include "force_volume.h"
#include "force_volume_yolk.h"
#include "equation_of_motion.h"
#include "dissolve.h"
#include "list_manipulation.h"
#include "T1_transformation.h"
#include "cell_division.h"
#include "cell_extrusion.h"
#include "LBFGS_helpers.h"
// //****************************************************************************
// //****************************************************************************