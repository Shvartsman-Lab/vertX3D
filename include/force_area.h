#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double f_SurfaceTension_force(const int i, int **_f){
    
    const unsigned int v1 = _f[i][1];
    const unsigned int v2 = _f[i][2];
    const unsigned int v3 = _f[i][3];

    //COORDINATES OF VERTICES
    //v3
    double v3x=v_pass[v3][1];
    double v3y=v_pass[v3][2];
    double v3z=v_pass[v3][3];
    
    //v1
    double v1x=v[v1][1];
    double v1y=v[v1][2];
    double v1z=v[v1][3];
    //x
    if(fabs(v1x-v3x)>0.5*perioXYZ[0])
        v1x += std::copysign(perioXYZ[0], v3x-v1x);
    //y
    if(fabs(v1y-v3y)>0.5*perioXYZ[1])
        v1y += std::copysign(perioXYZ[1], v3y-v1y);
    //z
    if(fabs(v1z-v3z)>0.5*perioXYZ[2])
        v1z += std::copysign(perioXYZ[2], v3z-v1z);
    
    //v2
    double v2x=v[v2][1];
    double v2y=v[v2][2];
    double v2z=v[v2][3];
    //x
    if(fabs(v2x-v3x)>0.5*perioXYZ[0])
        v2x += std::copysign(perioXYZ[0], v3x-v2x);
    //y
    if(fabs(v2y-v3y)>0.5*perioXYZ[1])
        v2y += std::copysign(perioXYZ[1], v3y-v2y);
    //z
    if(fabs(v2z-v3z)>0.5*perioXYZ[2])
        v2z += std::copysign(perioXYZ[2], v3z-v2z);
    
    
    //FACET AREA -- zapis z vektorskimi operacijami
    const double ax = v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y;
    const double ay = v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z;
    const double az = v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z);
    const double facet_area = 0.5*sqrt(pow(ax, 2.) + pow(ay, 2.) + pow(az, 2.));
    
    
    double grad_mult = (2e-3 * log(facet_area) + facet_area) / facet_area;
    //GRADIENT : te gradiente zapisi z vektorskimi operacijami!
    const double c0 = -0.25*f_T[i]*grad_mult/facet_area;
    // const double c0 = -0.25*f_T[i]/facet_area;
    const double c1 = (v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y);
    const double c2 = (v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z));
    const double c3 = (v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z);

    
    #pragma omp atomic
    v_F[v1][1] += c0*((-v2y + v3y)*c1 + (-v2z + v3z)*c2);
    #pragma omp atomic
    v_F[v1][2] += c0*(( v2x - v3x)*c1 + (-v2z + v3z)*c3);
    #pragma omp atomic
    v_F[v1][3] += c0*(( v2y - v3y)*c3 + ( v2x - v3x)*c2);
                  
    #pragma omp atomic
    v_F[v2][1] += c0*(( v1y - v3y)*c1 + ( v1z - v3z)*c2);
    #pragma omp atomic
    v_F[v2][2] += c0*((-v1x + v3x)*c1 + ( v1z - v3z)*c3);
    #pragma omp atomic
    v_F[v2][3] += c0*((-v1y + v3y)*c3 + (-v1x + v3x)*c2);

    
    return facet_area;
    
}
//****************************************************************************