#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//***************************VOLUMES******************************************
//****************************************************************************
//****************************************************************************
double dV(int i, int vert_ref_id, int **_f){//VOLUME OF FACET i
    
    //IDENTIFY VERTICES
    const int v1=_f[i][1];
    const int v2=_f[i][2];
    const int v3=_f[i][3];
    
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
    //vref
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    double vrefz = v[vert_ref_id][3];
    
    //PERIODIC BOUNDARY CONDITIONS
    double half_periox = 0.5*perioXYZ[0];
    double half_perioy = 0.5*perioXYZ[1];
    double half_perioz = 0.5*perioXYZ[2];
    //v1
    if(fabs(v1x-vrefx)>half_periox)
        v1x += std::copysign(perioXYZ[0], vrefx-v1x);
    if(fabs(v1y-vrefy)>half_perioy)
        v1y += std::copysign(perioXYZ[1], vrefy-v1y);
    if(fabs(v1z-vrefz)>half_perioz)
        v1z += std::copysign(perioXYZ[2], vrefz-v1z);
    //v2
    if(fabs(v2x-vrefx)>half_periox)
        v2x += std::copysign(perioXYZ[0], vrefx-v2x);
    if(fabs(v2y-vrefy)>half_perioy)
        v2y += std::copysign(perioXYZ[1], vrefy-v2y);
    if(fabs(v2z-vrefz)>half_perioz)
        v2z += std::copysign(perioXYZ[2], vrefz-v2z);
    //v3
    if(fabs(v3x-vrefx)>half_periox)
        v3x += std::copysign(perioXYZ[0], vrefx-v3x);
    if(fabs(v3y-vrefy)>half_perioy)
        v3y += std::copysign(perioXYZ[1], vrefy-v3y);
    if(fabs(v3z-vrefz)>half_perioz)
        v3z += std::copysign(perioXYZ[2], vrefz-v3z);

    //VOLUME
    return (-v1z*v2y*v3x + v1y*v2z*v3x + v1z*v2x*v3y - v1x*v2z*v3y - v1y*v2x*v3z + v1x*v2y*v3z)/6.0;
}
//****************************************************************************
double cell_volume(const int i, int **_f){
    
    int vert_ref_id=f[basal_facets[i][3]][1];
    double _cvol=0.0;
    
    #pragma omp simd reduction(+:_cvol)
    for(size_t j = 3; j <= 2+basal_edges[i][2]; ++j){
        //BASAL
        int fac_id=basal_facets[i][j];
        double sign = -1;
        _cvol += sign*dV(fac_id,vert_ref_id, _f);
        
        //APICAL
        fac_id = apical_facets[i][j];
        sign=1;
        _cvol += sign*dV(fac_id,vert_ref_id, _f);
        
        //LATERAL
        sign = std::copysign(1., basal_edges[i][j]);
        int edg_id = abs(basal_edges[i][j]);
        double _evol=0.0;
        //1
        fac_id=e_lateral1[edg_id];
        _evol += dV(fac_id,vert_ref_id, _f);
        //2
        fac_id=e_lateral2[edg_id];
        _evol += dV(fac_id,vert_ref_id, _f);
        //3
        fac_id=e_lateral3[edg_id];
        _evol += dV(fac_id,vert_ref_id, _f);
        //4
        fac_id=e_lateral4[edg_id];
        _evol += dV(fac_id,vert_ref_id, _f);
        
        _cvol += sign*_evol;
    }
    
    return _cvol;
}
//****************************************************************************
void f_VolCompressibility_force(int i, int cell_id, double sign, int vert_ref_id, double cell_vol, int **_f){

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
    //vref
    double vrefx = v[vert_ref_id][1];
    double vrefy = v[vert_ref_id][2];
    double vrefz = v[vert_ref_id][3];
    
    //PERIODIC BOUNDARY CONDITIONS
    //v1
    //x
    if(fabs(v1x-vrefx)>perioXYZ[0]/2.){
        if(v1x<vrefx) v1x+=perioXYZ[0];
        else if(v1x>vrefx) v1x-=perioXYZ[0];
    }
    //y
    if(fabs(v1y-vrefy)>perioXYZ[1]/2.){
        if(v1y<vrefy) v1y+=perioXYZ[1];
        else if(v1y>vrefy) v1y-=perioXYZ[1];
    }
    //z
    if(fabs(v1z-vrefz)>perioXYZ[2]/2.){
        if(v1z<vrefz) v1z+=perioXYZ[2];
        else if(v1z>vrefz) v1z-=perioXYZ[2];
    }
    
    //v2
    //x
    if(fabs(v2x-vrefx)>perioXYZ[0]/2.){
        if(v2x<vrefx) v2x+=perioXYZ[0];
        else if(v2x>vrefx) v2x-=perioXYZ[0];
    }
    //y
    if(fabs(v2y-vrefy)>perioXYZ[1]/2.){
        if(v2y<vrefy) v2y+=perioXYZ[1];
        else if(v2y>vrefy) v2y-=perioXYZ[1];
    }
    //z
    if(fabs(v2z-vrefz)>perioXYZ[2]/2.){
        if(v2z<vrefz) v2z+=perioXYZ[2];
        else if(v2z>vrefz) v2z-=perioXYZ[2];
    }
    
    //v3
    if(fabs(v3x-vrefx)>perioXYZ[0]/2.){
        if(v3x<vrefx) v3x+=perioXYZ[0];
        else if(v3x>vrefx) v3x-=perioXYZ[0];
    }
    //y
    if(fabs(v3y-vrefy)>perioXYZ[1]/2.){
        if(v3y<vrefy) v3y+=perioXYZ[1];
        else if(v3y>vrefy) v3y-=perioXYZ[1];
    }
    //z
    if(fabs(v3z-vrefz)>perioXYZ[2]/2.){
        if(v3z<vrefz) v3z+=perioXYZ[2];
        else if(v3z>vrefz) v3z-=perioXYZ[2];
    }
    
    double _deriv = -2.*sign*c_kV[cell_id]*(cell_vol-c_V0[cell_id]);
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
double c_VolCompressibility_force(const int cell_id, int **_f){
    
    int vert_ref_id=f[basal_facets[cell_id][3]][1];
    
    //CELL VOLUME
    double cell_vol=cell_volume(cell_id, _f);
    
    for(size_t j=3; j <= 2+basal_edges[cell_id][2]; j++){
        //BASAL
        f_VolCompressibility_force(basal_facets[cell_id][j],  cell_id, -1, vert_ref_id, cell_vol, _f);
        //APICAL
        f_VolCompressibility_force(apical_facets[cell_id][j], cell_id,  1, vert_ref_id, cell_vol, _f);
        
        //LATERAL
        int edg_id=abs(basal_edges[cell_id][j]);
        double sign = std::copysign(1,basal_edges[cell_id][j]);
        f_VolCompressibility_force(e_lateral1[edg_id], cell_id, sign, vert_ref_id, cell_vol, _f);
        f_VolCompressibility_force(e_lateral2[edg_id], cell_id, sign, vert_ref_id, cell_vol, _f);
        f_VolCompressibility_force(e_lateral3[edg_id], cell_id, sign, vert_ref_id, cell_vol, _f);
        f_VolCompressibility_force(e_lateral4[edg_id], cell_id, sign, vert_ref_id, cell_vol, _f);
    }
    return c_kV[cell_id]*pow((cell_vol-c_V0[cell_id]),2);
}
//****************************************************************************
