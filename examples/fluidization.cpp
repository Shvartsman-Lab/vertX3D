//****************GLOBAL PARAMETERS************************
double bet  =1.;    //basal tension
double alph =1.;    //apical tension
double kV   =100;   //reciprocal isothermal compressibility
double V0   =1.;    //preferred volume

//**************** SPECIFIC GLOBAL PARAMETERS************************
double kT1 = 50.; //Rate of Active T1 Transitions


//*******************INCLUDE*************************
#include "functions.h"
//***************************************************

void T1_merge(int i){
    if(basal_edges[e_cell1[i]][2]>3 && basal_edges[e_cell2[i]][2]>3 && v_cell4[e[i][1]]==0 && v_cell4[e[i][2]]==0 ){
        //printf("T1 on edge %d\n", i);
        int vertid=merge_vertices(i);
        v_clock[vertid]=0;
    }
    /*else if (basal_edges[e_cell1[i]][2]==3 || basal_edges[e_cell2[i]][2]==3){
        if (basal_edges[e_cell1[i]][2]==3) printf("T1 on edge %d ... can't perform (triangular cell %d)\n",i, e_cell1[i]);
        if (basal_edges[e_cell2[i]][2]==3) printf("T1 on edge %d ... can't perform (triangular cell %d)\n",i, e_cell2[i]);
    }
    else if (v_cell4[e[i][1]]!=0 || v_cell4[e[i][2]]!=0){
        if (v_cell4[e[i][1]]!=0) printf("T1 on edge %d ... can't perform (4-way vertex %d)\n", i, e[i][1]);
        if (v_cell4[e[i][2]]!=0) printf("T1 on edge %d ... can't perform (4-way vertex %d)\n", i, e[i][2]);
    }*/
}
//****************************************************************************
void val_red(double vclock, double finLength){
    
    for(size_t i = 1; i<=Nv; ++i){
        if( v[i][0]>0.5 && v_cell4[i]!=0 ){
            if(v_clock[i]>vclock){
                //printf("valence reduction on vertex %d\n", (int)i);
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
                //std::clog << "Spontaneous ";
                T1_merge(i);
            }
            //ACTIVE
            else if( rnd()<_kT1*h/(1.*Ne) ){
                //std::clog << "Active ";
                T1_merge(i);
            }
            else e_length[i]=eLen;
        }
    }
    
}
//****************************************************************************
/*double ceHeight(int i){
    
    double ddx,ddy,ddz;
    double *dxdydz = new double[3];
    dxdydz[0]=0; dxdydz[1]=0; dxdydz[2]=0;
    
    int v1=c_cent_apical[i], v2=c_cent_basal[i];
    torus_pass_dx_dy_dz(dxdydz,v1,v2);
    ddx=v_pass[v2][1]-(v_pass[v1][1]+dxdydz[0]);
    ddy=v_pass[v2][2]-(v_pass[v1][2]+dxdydz[1]);
    ddz=v_pass[v2][3]-(v_pass[v1][3]+dxdydz[2]);
    double hght=sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
    
    delete []dxdydz;
    
    return hght;
}*/
//****************************************************************************
void outDeformationField(int intTime){
    
    char fileName[50];
    snprintf(fileName, sizeof(char) * 50, "./output/defField_%d.txt", intTime);
    FILE *file1; file1 = fopen(fileName, "wt");
    
    double avgz=0;
    for(int i=1; i<=Nc; i++) if(basal_edges[i][1]!=0) avgz+=v_pass[i][3]/(1.*Nc);
    
    int vertID;
    double cellHeight;
    for(int i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0) {
            vertID=c_cent_apical[i];
            fprintf(file1, "%g  %g  %g  %g\n", v_pass[vertID][1], v_pass[vertID][2], v_pass[vertID][3], avgz);
        }
    }
    
    fclose(file1);
}
//****************************************************************************
void run(){
    
    set_initial_fromFile("./initial/crumpled340.vt3d");
    //*************************
    
    //*********************************************************************************
    //*********************************************************************************
    //*********************************************************************************
    //SIMULATION
    //*************************
    //DYNAMICS
    //expand_box(-0.15,-0.15,0);
    FILE *filee;
    filee = fopen("./output/energy.txt", "wt"); fclose(filee);
    
    double tcount=100, tcount2=100;
    int timeINT=0;
    
    h=0.001;
    double hess_guess = 1.0;
    int do_implicit = 0;
    bool yolk_present = false;
    while(Time<tmax){
        
        if(tcount>10){
            timeINT++;
            outDeformationField(timeINT);
            // outStructure(timeINT);
            write_vtk(timeINT);
            tcount=0;
        }
        
        //EQUATION OF MOTION
        if (eqOfMotion(hess_guess, do_implicit, yolk_present) != 0) {
            printf("Exiting time loop\n");
            break;
        }
        
        if(tcount2>1){
            filee = fopen("./output/energy.txt", "at");
            fprintf(filee, "%g  %g  %g\n", Time, wA, wV);
            tcount2=0;
            fclose(filee);
        }
        
        //TOPOLOGICAL TRANSFORMATIONS
        val_red(20*h,0.001);
        T1_spont_act(0.1,300-300*Time/tmax);
        // T1_spont_act(0.1,kT1);
        
        tcount+=h;
        tcount2+=h;
    }
    printf("%d  %d  %d  %d\n", (int)Nv, (int)Ne, (int)Nc, (int)Nv_pass);
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
