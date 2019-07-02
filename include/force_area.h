#pragma once

#include "functions.h"

//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
double f_SurfaceTension_force(const int i, int **_f, double implicit_eps){
    
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
    

    //FACET AREA
    double ax = v2y*v1z - v1y*v2z + v1y*v3z - v3y*v1z - v2y*v3z + v3y*v2z;
    double ay = v1x*v2z - v2x*v1z - v1x*v3z + v3x*v1z + v2x*v3z - v3x*v2z;
    double az = v2x*v1y - v1x*v2y + v1x*v3y - v3x*v1y - v2x*v3y + v3x*v2y;
    double facet_area = 0.5 * sqrt(ax*ax + ay*ay + az*az);

    // double neohookean = (2e-3 * log(facet_area) + facet_area) / facet_area;
    double grad_mult = 0.25 * f_T[i] / (facet_area + implicit_eps);
    // double grad_mult = 0.25 * f_T[i] / facet_area;
    int S_m = 4;
    if (f_type[i] != 3) { //Check if lateral
        S_m = basal_edges[f_cell[i]][2];
    }

    double a1 = az*(v2y - v3y) - ay*(v2z - v3z);
    double a2 = ax*(v2z - v3z) - az*(v2x - v3x);
    double a3 = ay*(v2x - v3x) - ax*(v2y - v3y);

    double b1 = az*(v1y - v3y) - ay*(v1z - v3z);
    double b2 = ax*(v1z - v3z) - az*(v1x - v3x);
    double b3 = ay*(v1x - v3x) - ax*(v1y - v3y);

    #pragma omp atomic
    v_F[v1][1] += grad_mult * a1;
    #pragma omp atomic
    v_F[v1][2] += grad_mult * a2;
    #pragma omp atomic
    v_F[v1][3] += grad_mult * a3;

    #pragma omp atomic
    v_F[v2][1] += grad_mult * (-b1);
    #pragma omp atomic
    v_F[v2][2] += grad_mult * (-b2);
    #pragma omp atomic
    v_F[v2][3] += grad_mult * (-b3);

    double c1 = grad_mult * (-(1.0 / S_m) * a1 + (1.0 / S_m) * b1);
    double c2 = grad_mult * (-(1.0 / S_m) * a2 + (1.0 / S_m) * b2);
    double c3 = grad_mult * (-(1.0 / S_m) * a3 + (1.0 / S_m) * b3);

    int v_id;
    double d_ip, d_iq;
    if (f_type[i] != 3) {
        //IF BASAL OR APICAL
        for (size_t j = 3; j <= 2 + S_m; j++) {
            if (f_type[i] == 1) {
                v_id = basal_vertices[f_cell[i]][j];
            } else {
                v_id = v_partner[basal_vertices[f_cell[i]][j]];
            }
            
            #pragma omp atomic
            v_F[v_id][1] += c1;
            #pragma omp atomic
            v_F[v_id][2] += c2;
            #pragma omp atomic
            v_F[v_id][3] += c3;
        }
    } else {
        //IF LATERAL
        //lateral1
        // v_id = _f[e_lateral1[edg_id]][1];
        v_id = e[f_edge[i]][1];
        #pragma omp atomic
        v_F[v_id][1] += c1;
        #pragma omp atomic
        v_F[v_id][2] += c2;
        #pragma omp atomic
        v_F[v_id][3] += c3;

        //lateral2
        // v_id = _f[e_lateral2[edg_id]][1];
        v_id = e[f_edge[i]][2];
        #pragma omp atomic
        v_F[v_id][1] += c1;
        #pragma omp atomic
        v_F[v_id][2] += c2;
        #pragma omp atomic
        v_F[v_id][3] += c3;

        //lateral3
        // v_id = _f[e_lateral3[edg_id]][1];
        v_id = v_partner[e[f_edge[i]][1]];
        #pragma omp atomic
        v_F[v_id][1] += c1;
        #pragma omp atomic
        v_F[v_id][2] += c2;
        #pragma omp atomic
        v_F[v_id][3] += c3;

        //lateral4
        // v_id = _f[e_lateral4[edg_id]][1];
        v_id = v_partner[e[f_edge[i]][2]];
        #pragma omp atomic
        v_F[v_id][1] += c1;
        #pragma omp atomic
        v_F[v_id][2] += c2;
        #pragma omp atomic
        v_F[v_id][3] += c3;
    }

    return facet_area;
}
//****************************************************************************