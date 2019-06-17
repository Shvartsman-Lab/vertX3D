//****************GLOBAL PARAMETERS************************
double bet  =1.;    //basal tension
double alph =1.;    //apical tension
double kV   =100;   //reciprocal isothermal compressibility
double V0   =1.;    //preferred volume

//*******************INCLUDE*************************
#include "functions.h"
//***************************************************

double dA(int i){    
    const unsigned int v1 = f[i][1];
    const unsigned int v2 = f[i][2];
    const unsigned int v3 = f[i][3];
    
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
    const double ax = v1y*v2x - v1x*v2y - v1y*v3x + v2y*v3x + v1x*v3y - v2x*v3y;
    const double ay = v1z*v2y - v1y*v2z - v1z*v3y + v2z*v3y + v1y*v3z - v2y*v3z;
    const double az = v1z*(v2x - v3x) + v2z*v3x - v2x*v3z + v1x*(-v2z + v3z);
    
    return 0.5*sqrt(pow(ax, 2.) + pow(ay, 2.) + pow(az, 2.));
    
}
//****************************************************************************
double Amid(int i){
    
    double Asum_a, Asum_b;
    if(basal_edges[i][1]>0){
        Asum_a=0;
        Asum_b=0;
        for(size_t j = 3; j<=2+basal_edges[i][2]; ++j){
            Asum_a+=dA(apical_facets[i][j]);
            Asum_b+=dA(basal_facets[i][j]);
        }
    }
    
    return 0.5*(Asum_a+Asum_b);
}
//****************************************************************************
void project_on_constraints(double z0){
    
    for(size_t i=1; i<=Nv; i++){
        if(v[i][0]>0.5 && v_type[i]==1){
            v[i][3]=z0;
        }
    }

}
//****************************************************************************
void T1_merge(int i){
    if(basal_edges[e_cell1[i]][2]>3 && basal_edges[e_cell2[i]][2]>3 && v_cell4[e[i][1]]==0 && v_cell4[e[i][2]]==0 ){
        printf("T1 on edge %d\n", i);
        int vertid=merge_vertices(i);
        v_clock[vertid]=0;
    }
    else if (basal_edges[e_cell1[i]][2]==3 || basal_edges[e_cell2[i]][2]==3){
        if (basal_edges[e_cell1[i]][2]==3) printf("T1 on edge %d ... can't perform (triangular cell %d)\n",i, e_cell1[i]);
        if (basal_edges[e_cell2[i]][2]==3) printf("T1 on edge %d ... can't perform (triangular cell %d)\n",i, e_cell2[i]);
    }
    else if (v_cell4[e[i][1]]!=0 || v_cell4[e[i][2]]!=0){
        if (v_cell4[e[i][1]]!=0) printf("T1 on edge %d ... can't perform (4-way vertex %d)\n", i, e[i][1]);
        if (v_cell4[e[i][2]]!=0) printf("T1 on edge %d ... can't perform (4-way vertex %d)\n", i, e[i][2]);
    }
}
//****************************************************************************
void val_red(double vclock, double finLength){
    
    for(size_t i = 1; i<=Nv; ++i){
        if( v[i][0]>0.5 && v_cell4[i]!=0 ){
            if(v_clock[i]>vclock){
                valence_reduction(i,finLength);
            }
            else v_clock[i]+=h;
        }
    }
}
//****************************************************************************
void T1_spont_act(double _threshold_length, double _kT1){
    
    for(size_t i = 1; i<=Ne; ++i){
        if( e[i][0]!=0 ){
            
            double eLen=edge_length(i);
            
            //SPONTANEOUS
            if( eLen<_threshold_length && eLen<e_length[i]){
                std::clog << "Spontaneous ";
                T1_merge(i);
            }
            //ACTIVE
            else if( rnd()<_kT1*h/(1.*Ne) ){
                std::clog << "Active ";
                T1_merge(i);
            }
            else e_length[i]=eLen;
        }
    }
    
}
//****************************************************************************
int division_extrusion(double kd, double ke, double Temp){
    
    double dw, Am, et;
    int cell_balance=0;
    
    for(size_t i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            
            Am=Amid(i);
            et=basal_edges[i][2]*tan(M_PI/(1.*basal_edges[i][2]))/4.;
            dw=(alph+bet)*Am+2.*sqrt(et/Am)-3.*pow(et*(alph+bet),1./3.);
            
            //EXTRUSION
            if ( Am<pow(et,1./3.)*pow(alph+bet,-2./3.) && rnd()<ke*(1.-exp(-dw/Temp)) ) cell_balance-=cell_extrusion(i);
            //DIVISION
            else if ( Am>pow(et,1./3.)*pow(alph+bet,-2./3.) && rnd()<kd*(1.-exp(-dw/Temp)) ) cell_balance+=cell_division(i);
        }
    }
    return cell_balance;
}
//****************************************************************************
void run(){    
    set_initial_fromFile("./initial/stretched.vt3d");
    char filename[50]; snprintf(filename, sizeof(char) * 50, "./output/out_NrOfCells.txt");
    FILE *fNrCells; fNrCells = fopen(filename, "wt");
    h=0.001;
    int nr_of_cells=Nc;

    double hess_guess = 1.0;
    int do_implicit = 0;
    while(Time<tmax){
        
        //EQUATION OF MOTION
        if (eqOfMotion(hess_guess, do_implicit) != 0) {
            printf("Exiting time loop\n");
            break;
        }
        project_on_constraints(5.);
        
        //VALENCE REDUCTION ON 4-WAY VERTICES
        val_red(20*h,0.001);
        
        //T1 TRANSFORMATION
        T1_spont_act(0.15,0);
        
        //DIVISION/EXTRUSION
        nr_of_cells+=division_extrusion(0.0002,0.0002,0.03);
        fprintf(fNrCells,"%g  %d\n", Time, nr_of_cells);
    }
    // outStructure(0);
    write_vtk(0);
    
    fclose(fNrCells);
    
}

//********************MAIN********************************
int main(int argc, char *argv[]){
    read_code_arguments(argc,argv);
    allocate();
    run();
    out_Vertissue3D("./output/out_Final.vt3d");
    deallocate();
    return 0;
}
//********************************************************
