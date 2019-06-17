#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***********************INITIAL STRUCTURE************************************
//****************************************************************************
//****************************************************************************
void set_initial_regHex(size_t _Nx){
    
    //INTERNAL VARIABLES
    double *v_stitching_edge;//only useful for building regular hexagonal network
    v_stitching_edge = new double[array_max];
    int Ny=_Nx/2;
    double l0=pow(((1+1)*(1+1))/(0.5*sqrt(3.)),0.333333);
    double dd=(2./sqrt(3.))*pow(0.5*sqrt(3.)*(1+1),-0.333333);
    double z0=5., ddx=0.08;
    
    //GLOBAL VARIABLES
    perioXYZ[0] = _Nx*sqrt(3)*dd/2.;
    perioXYZ[1] = Ny*(dd+dd/2);
    perioXYZ[2] = 10;
    
    for(size_t i=1; i<=Ny; i++){
        for (size_t j=1; j<=_Nx; j++){
            double x0=(j-1)*sqrt(3)*dd/2+ddx;
            double y0=(3*dd/2)*(i-1)+ddx;
            make_vertex(
                        x0,
                        y0-dd/2+dd/2,
                        z0
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0-dd/4+dd/2,
                        z0
                        );
            make_vertex(
                        x0+sqrt(3)*dd/4,
                        y0+dd/4+dd/2,
                        z0
                        );
            make_vertex(
                        x0,
                        y0+dd/2+dd/2,
                        z0
                        );
        }
        
        make_edge(i*_Nx*4-2, (i-1)*_Nx*4+1);
        make_edge((i-1)*_Nx*4+1, (i-1)*_Nx*4+2);
        make_edge((i-1)*_Nx*4+2, (i-1)*_Nx*4+3);
        make_edge((i-1)*_Nx*4+3, (i-1)*_Nx*4+4);
        make_edge((i-1)*_Nx*4+4, i*_Nx*4-1);
        
        for (size_t j=2; j<=_Nx; j++){
            make_edge((i-1)*_Nx*4+(j-2)*4+2, (i-1)*_Nx*4+(j-2)*4+5);
            make_edge((i-1)*_Nx*4+(j-2)*4+5, (i-1)*_Nx*4+(j-2)*4+6);
            make_edge((i-1)*_Nx*4+(j-2)*4+6, (i-1)*_Nx*4+(j-2)*4+7);
            make_edge((i-1)*_Nx*4+(j-2)*4+7, (i-1)*_Nx*4+(j-2)*4+8);
            make_edge((i-1)*_Nx*4+(j-2)*4+8, (i-1)*_Nx*4+(j-2)*4+3);
        }
    }
    
    for (size_t i=1; i<=Ny-1; i++){
        for (size_t j=1; j<=_Nx; j++){
            int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, _Nx*4*i+1+(j-1)*4);
            v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
        }
    }
    
    int i = Ny;
    for (size_t j=1; j<=_Nx; j++){
        int edgid=make_edge(_Nx*4*(i-1)+4+(j-1)*4, 1+(j-1)*4);
        v_stitching_edge[_Nx*4*(i-1)+4+(j-1)*4]=edgid;
    }
    
    
    
    //*************************************
    //*************************************
    //*************************************
    //BASAL EDGES
    //*************************************x   
    //LIHE VRSTICE
    for (size_t i=1; i<=Ny; i++){
        
        make_cell(
        6,
        (i-1)*_Nx*5+1,
        (i-1)*_Nx*5+2,
        (i-1)*_Nx*5+3,
        (i-1)*_Nx*5+4,
        (i-1)*_Nx*5+5,
        -((i-1)*_Nx*5+(_Nx*5-2)),
        0,0,0,0,0,0
        );
        
        for (size_t j=6; j<=_Nx*5-4; j+=5){
            make_cell(
            6,
            (i-1)*_Nx*5+j,
            (i-1)*_Nx*5+j+1,
            (i-1)*_Nx*5+j+2,
            (i-1)*_Nx*5+j+3,
            (i-1)*_Nx*5+j+4,
            -((i-1)*_Nx*5+j-3),
            0,0,0,0,0,0
            );
        }
    }
    
    //SODE VRSTICE
    for (size_t i=1; i<=Ny-1; i++){
        
        for (size_t j=1; j<=_Nx-1; j++){
            make_cell(
            6,
            -(4+(j-1)*5+(i-1)*_Nx*5),
            -(10+(j-1)*5+(i-1)*_Nx*5),
            v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]],
            -(5*_Nx+6+(j-1)*5+(i-1)*_Nx*5),
            -(5*_Nx+6+(j-1)*5-4+(i-1)*_Nx*5),
            -(v_stitching_edge[e[(10+(j-1)*5)+(i-1)*_Nx*5][1]]-1),
            0,0,0,0,0,0
            );
        }
        
        make_cell(
        6,
        -(4+(_Nx-1)*5+(i-1)*_Nx*5),
        -(5+(i-1)*_Nx*5),
        v_stitching_edge[e[10+(i-1)*_Nx*5][1]]-1,
        -(5*_Nx+6-4-1+(i-1)*_Nx*5),
        -(5*_Nx+6+(_Nx-1)*5-4+(i-1)*_Nx*5),
        -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(i-1)*_Nx*5)][1]]),
        0,0,0,0,0,0
        );
    }
    
    //ZADNJA VRSTICA
    for (size_t j=1; j<=_Nx-1; j++){
        
        make_cell(
        6,
        -(4+(j-1)*5+(Ny-1)*_Nx*5),
        -(10+(j-1)*5+(Ny-1)*_Nx*5),
        v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]],
        -(6+(j-1)*5),
        -(2+(j-1)*5),
        -(v_stitching_edge[e[(10+(j-1)*5)+(Ny-1)*_Nx*5][1]]-1),
        0,0,0,0,0,0
        );
    }
    
    make_cell(
    6,
    -(4+(_Nx-1)*5+(Ny-1)*_Nx*5),
    -(5+(Ny-1)*_Nx*5),
    v_stitching_edge[e[10+(Ny-1)*_Nx*5][1]]-1,
    -1,
    -(_Nx*5-3),
    -(v_stitching_edge[e[(10+(_Nx-1-1)*5+(Ny-1)*_Nx*5)][1]]),
    0,0,0,0,0,0
    );
    
    
    //*****************************
    //*****************************
    //*****************************
    for(size_t i=1; i<=Nv; i++){
        v_normal_vector[i][1]=0;
        v_normal_vector[i][2]=0;
        v_normal_vector[i][3]=1;
        v_height[i]=l0;
    }
    
    
    delete [] v_stitching_edge;
    
    //CREATE 3D TISSUE
    make_basal_side_ALL();
    copy_vertex_to_apical_ALL();
    make_apical_side_ALL();
    make_lateral_side_ALL();
    
    std::clog << "Initialization completed." << std::endl;
}
//****************************************************************************
void set_initial_fromFile(char *fileName){
    
    FILE *file1;
    file1 = fopen(fileName, "rt");
    
    //V,E,C
    int nrV=0, nrE=0, nrC=0;
    fscanf(file1, "%d  %d  %d\n", &nrV, &nrE, &nrC);
    
    //perioXYZ
    fscanf(file1, "%lf  %lf  %lf\n", &perioXYZ[0], &perioXYZ[1], &perioXYZ[2]);
    
    //printf("%d  %d  %d\n", nrV, nrE, nrC);
    //printf("%g  %g  %g\n", perioXYZ[0], perioXYZ[1], perioXYZ[2]);
    
    //VERTICES
    double xx, yy, zz, nx, ny, nz, hh;
    for(size_t i=1; i<=nrV; i++){
        fscanf(file1, "%lf  %lf  %lf  %lf  %lf  %lf  %lf\n", &xx, &yy, &zz, &nx, &ny, &nz, &hh);
        make_vertex(xx,yy,zz);
        v_normal_vector[i][1]=nx;
        v_normal_vector[i][2]=ny;
        v_normal_vector[i][3]=nz;
        v_height[i]=hh;
    }
    
    //EDGES
    int v1, v2;
    for(size_t i=1; i<=nrE; i++){
        fscanf(file1, "%d  %d\n", &v1, &v2);
        make_edge(v1,v2);
    }
    
    //CELLS
    int b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13;
    for(size_t i=1; i<=nrC; i++){
        fscanf(file1, "%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n", &b1, &b2, &b3, &b4, &b5, &b6, &b7, &b8, &b9, &b10, &b11, &b12, &b13);
        make_cell(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13);
    }
    fclose(file1);
    
    //CREATE 3D TISSUE
    make_basal_side_ALL();
    copy_vertex_to_apical_ALL();
    make_apical_side_ALL();
    make_lateral_side_ALL();
    
    std::clog << "Initialization completed." << std::endl;
}
//****************************************************************************
