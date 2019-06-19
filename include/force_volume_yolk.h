#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************VOLUMES******************************************
//****************************************************************************
//****************************************************************************

double dVY(int i, int **_f){//VOLUME OF FACET i
    
    //IDENTIFY VERTICES
    const int v1 =_f[i][1];
    const int v2 =_f[i][2];
    const int v3 =_f[i][3];
    
    //COORDINATES OF VERTICES
    //v1
    double v1x = v[v1][1];
    double v1y = v[v1][2];
    double v1z = v[v1][3];
    //v2
    double v2x = v[v2][1];
    double v2y = v[v2][2];
    double v2z = v[v2][3];
    //v3
    double v3x = v_pass[v3][1];
    double v3y = v_pass[v3][2];
    double v3z = v_pass[v3][3];
    
    //VOLUME
    return (-v1z*v2y*v3x + v1y*v2z*v3x + v1z*v2x*v3y - v1x*v2z*v3y - v1y*v2x*v3z + v1x*v2y*v3z)/6.0;
}
//****************************************************************************
double yolk_volume(int **_f){
    
    double _Yvol =0.0;
    int sign=1, fac_id;
    
    for(size_t i=1; i<=Nc;i++){
        if(basal_edges[i][1]!=0){
            for(size_t j = 3; j <= 2+basal_edges[i][2]; ++j){
                fac_id = apical_facets[i][j];
                _Yvol += sign*dVY(fac_id, _f);
            }
        }
    }
    
    return _Yvol;
}
//****************************************************************************
void f_YolkCompressibility_force(int i, double sign, double Yolk_vol, int **_f){
    
    //IDENTIFY VERTICES
    const int v1 = _f[i][1];
    const int v2 = _f[i][2];
    const int v3 = _f[i][3];
    
    //COORDINATES OF VERTICES
    //v1
    double v1x = v[v1][1];
    double v1y = v[v1][2];
    double v1z = v[v1][3];
    //v2
    double v2x = v[v2][1];
    double v2y = v[v2][2];
    double v2z = v[v2][3];
    //v3
    double v3x = v_pass[v3][1];
    double v3y = v_pass[v3][2];
    double v3z = v_pass[v3][3];
    
    double _deriv = -2.*sign*c_kY*(Yolk_vol-VY0);
    //GRADIENT
    //v1
    #pragma omp atomic
    v_F[v1][1] += _deriv*(-(v2z*v3y) + v2y*v3z)/6.;
    #pragma omp atomic
    v_F[v1][2] += _deriv*(  v2z*v3x  - v2x*v3z)/6.;
    #pragma omp atomic
    v_F[v1][3] += _deriv*(-(v2y*v3x) + v2x*v3y)/6.;
    //v2
    #pragma omp atomic
    v_F[v2][1] += _deriv*(  v1z*v3y  - v1y*v3z)/6.;
    #pragma omp atomic
    v_F[v2][2] += _deriv*(-(v1z*v3x) + v1x*v3z)/6.;
    #pragma omp atomic
    v_F[v2][3] += _deriv*(  v1y*v3x  - v1x*v3y)/6.;
    
}
//****************************************************************************
double YolkCompressibility_force(int **_f){
    
    //CELL VOLUME
    double _Yvol = yolk_volume(_f);
    int sign=1;
    
    #pragma omp parallel for schedule(guided) shared(v_F, v)
    for(size_t i=1; i<=Nc;i++){
        if(basal_edges[i][1]!=0){
            for(size_t j =3; j <= 2+basal_edges[i][2]; j++){
                int fac_id=apical_facets[i][j];
                f_YolkCompressibility_force(fac_id, sign, _Yvol, _f);
            }
        }
    }
    return c_kY*pow(_Yvol-VY0,2);
}
//****************************************************************************