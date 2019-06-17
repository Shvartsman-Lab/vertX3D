#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//******************************ALLOCATE**************************************
//****************************************************************************
//****************************************************************************
int read_code_arguments(int argc, char *argv[]){
    
    //checks if there are enough arguments specified
    if(argc != 5){
        std::cerr << "Error! Wrong number of input elements!" << std::endl;
        std::cerr << "Run as: prog.out array_max tmax seed omp_num_threads\nTerminating" << std::endl;
        return -1;
    }
    
    //array_size
    array_max = std::stol(argv[1]);
    std::clog << "array_max: " << array_max << std::endl;
    
    //tmax
    tmax = std::stod(argv[2]);
    std::clog << "tmax: " << tmax << std::endl;
    
    //seed
    seed = std::stol(argv[3]);
    std::clog << "seed: " << seed << std::endl;
    srand(seed);
    
    //number of threads
    size_t num_threads=std::stoul(argv[4]);
    std::clog << "NumThreads: " << num_threads << std::endl;
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    
    return 1;
}
//****************************************************************************
void reset_arrays(){
    
    #pragma omp parallel for
    for(size_t i = 1; i<=array_max; i++){
        
        //***********************
        //GEOMETRIC ELEMENTS*****
        //***********************
        for(size_t j = 0; j < 4; j++){
            v[i][j]=0;
            v_pass[i][j]=0;
            f[i][j]=0;
        }
        for(size_t j = 0; j < 3; j++){
            e[i][j]=0;
        }
        //basal_edges
        basal_edges[i][1]=0;
        basal_edges[i][2]=0;
        //basal_vertices
        basal_vertices[i][1]=0;
        basal_vertices[i][2]=0;
        //basal_facets
        basal_facets[i][1]=0;
        basal_facets[i][2]=0;
        //apical_facets
        apical_facets[i][1]=0;
        apical_facets[i][2]=0;
        //***********************
        
        
        //***********************
        //ATTRIBUTES*************
        //***********************
        //vertices
        for(size_t j = 1; j < 4; j++){
            v_F[i][j]=0;
            v_normal_vector[i][j]=0;
        }
        v_height[i]=0;
        v_type[i]=0;
        v_partner[i]=0;
        v_cell1[i]=0;
        v_cell2[i]=0;
        v_cell3[i]=0;
        v_cell4[i]=0;
        v_T1dir[i]=0;
        v_vertT1[i]=0;
        v_edgeT1[i]=0;
        v_clock[i]=0;

        //passive vertices
        v_pass_type[i]=0;
        v_pass_cell[i]=0;
        v_pass_edge[i]=0;

        //edges
        e_v_pass[i]=0;
        e_lateral1[i]=0;
        e_lateral2[i]=0;
        e_lateral3[i]=0;
        e_lateral4[i]=0;
        e_cell1[i]=0;
        e_cell2[i]=0;
        e_length[i]=0;

        //facets
        f_type[i]=0;
        f_cell[i]=0;
        f_edge[i]=0;
        f_T[i]=0;

        //cells
        c_cent_basal[i]=0;
        c_cent_apical[i]=0;
        c_V0[i]=0;
        c_kV[i]=0;
        c_alpha[i]=0;
        c_beta[i]=0;
    }
}
//****************************************************************************
void allocate(){// allocate space
    
    size_t arrayMax=array_max+1;
    //std::clog << "Max Array Size: " << arrayMax << std::endl;
    
    //GEOMETRIC ELEMENTS
    Nv      = 0;
    Nv_pass = 0;
    Ne      = 0; 
    Nf      = 0; 
    Nc      = 0;
    v        = new double*[arrayMax];
    v_pass   = new double*[arrayMax];
    e = new int*[arrayMax];
    f = new int*[arrayMax];
    basal_edges    = new int*[arrayMax];
    basal_vertices = new int*[arrayMax];
    basal_facets   = new int*[arrayMax];
    apical_facets  = new int*[arrayMax];
    for(size_t i=0; i < arrayMax; i++){
        v[i]        = new double[4];
        v_pass[i]   = new double[4];
        e[i]   = new int[3];
        f[i]   = new int[4];
        basal_edges[i]    = new int[15];
        basal_vertices[i] = new int[15];
        basal_facets[i]   = new int[15];
        apical_facets[i]  = new int[15];
    }
    
    //ATTRIBUTES
    //vertices
    v_F      = new double*[arrayMax];
    v_normal_vector = new double*[arrayMax];
    for(size_t i=0; i < arrayMax; i++){
        v_F[i]      = new double[4];
        v_normal_vector[i] = new double[4];
    }
    v_height = new double[arrayMax];
    v_type     = new int[arrayMax];
    v_partner  = new int[arrayMax];
    v_cell1    = new int[arrayMax];
    v_cell2    = new int[arrayMax];
    v_cell3    = new int[arrayMax];
    v_cell4    = new int[arrayMax];
    v_T1dir    = new int[arrayMax];
    v_vertT1   = new int[arrayMax];
    v_edgeT1   = new int[arrayMax];
    v_clock   = new double[arrayMax];
    
    //passive vertices
    v_pass_type = new int[arrayMax];
    v_pass_cell = new int[arrayMax];
    v_pass_edge = new int[arrayMax];

    //edges
    e_v_pass   = new int[arrayMax];
    e_lateral1 = new int[arrayMax];
    e_lateral2 = new int[arrayMax];
    e_lateral3 = new int[arrayMax];
    e_lateral4 = new int[arrayMax];
    e_cell1    = new int[arrayMax];
    e_cell2    = new int[arrayMax];
    e_length = new double[arrayMax];

    //facets
    f_type = new int[arrayMax];
    f_cell = new int[arrayMax];
    f_edge = new int[arrayMax];
    f_T = new double[arrayMax];

    //cells
    c_cent_basal = new int[arrayMax];
    c_cent_apical = new int[arrayMax];
    c_V0 = new double[arrayMax];
    c_kV = new double[arrayMax];
    c_alpha = new double[arrayMax];
    c_beta = new double[arrayMax];
    
    //MISC
    perioXYZ = new double[3];
    h              = 0;
    Time           = 0;
    max_move       = 0;
    wA             = 0;
    wV             = 0;
    A_tot          = 0;
    
    //RESET ARRAYS
    reset_arrays();

    return ;
}
//****************************************************************************
void deallocate(){// Deallocate space
    
    size_t arrayMax=array_max+1;

    //GEOMETRIC ELEMENTS
    for(size_t i=0; i < arrayMax; i++){
        delete [] v[i]      ;
        delete [] v_pass[i];
        delete [] e[i];
        delete [] f[i];
        delete [] basal_edges[i]   ;
        delete [] basal_vertices[i];
        delete [] basal_facets[i]  ;
        delete [] apical_facets[i] ;
    }
    delete [] v;
    delete [] v_pass;
    delete [] e;
    delete [] f;
    delete [] basal_edges   ;
    delete [] basal_vertices;
    delete [] basal_facets  ;
    delete [] apical_facets ;

    //ATTRIBUTES
    //vertices
    for(size_t i=0; i < arrayMax; i++){
        delete [] v_F[i]    ;
        delete [] v_normal_vector[i];
    }
    delete [] v_F    ;
    delete [] v_normal_vector;
    delete [] v_height;
    delete [] v_type    ;
    delete [] v_partner ;
    delete [] v_cell1   ;
    delete [] v_cell2   ;
    delete [] v_cell3   ;
    delete [] v_cell4   ;
    delete [] v_T1dir   ;
    delete [] v_vertT1  ;
    delete [] v_edgeT1  ;
    delete [] v_clock  ;

    //passive vertices
    delete [] v_pass_type;
    delete [] v_pass_cell;
    delete [] v_pass_edge;
    
    //edges
    delete [] e_v_pass  ;
    delete [] e_lateral1;
    delete [] e_lateral2;
    delete [] e_lateral3;
    delete [] e_lateral4;
    delete [] e_cell1   ;
    delete [] e_cell2   ;
    delete [] e_length;

    
    //facets
    delete [] f_type;
    delete [] f_cell;
    delete [] f_edge;
    delete [] f_T;
    
    //cells
    delete [] c_cent_basal;
    delete [] c_cent_apical;
    delete [] c_V0;
    delete [] c_kV;
    delete [] c_alpha;
    delete [] c_beta;
    
    //MISC
    delete [] perioXYZ;

    return ;
}
//****************************************************************************