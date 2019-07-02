//****************GLOBAL PARAMETERS************************
double bet  =1.;    //basal tension
double alph =1.;    //apical tension
double kV   =100;   //reciprocal isothermal compressibility
double V0   =1.;    //preferred volume

//*******************INCLUDE*************************
#include "functions.h"
//***************************************************

void project_on_constraints(double r0){
    for(size_t i=1; i<=Nv; i++) if(v[i][0]>0.5 && v_type[i]==1) v[i][3]=r0;
}
//****************************************************************************
void spont_T1(double _threshold_length){
    
    for(size_t i = 1; i<=Ne; ++i){
        
        if( e[i][0]!=0 ){
            
            double eLen=edge_length(i);
            
            if( eLen<_threshold_length && eLen<e_length[i]){
                std::clog << "Spontaneous ";
                T1(i,2*_threshold_length);
            }
            
            else e_length[i]=eLen;
            
        }
    }
    
}
//****************************************************************************
int divide_cell(double k){
    int cell_balance=0;
    for(size_t i=1; i<=Nc; i++) if(basal_edges[i][1]!=0) if ( rnd()<k ) cell_balance+=cell_division(i);
    return cell_balance;
}
//****************************************************************************
void countPolygons(){
    
    FILE *file1; file1 = fopen("./output/polygonality.txt", "at");
    
    int *polygonality = new int[30];
    for(int i=3; i<=12; i++) polygonality[i]=0;
    
    int nr_of_cells=0;
    for(size_t i=1; i<=Nc; i++){
        if(basal_edges[i][1]!=0){
            polygonality[basal_edges[i][2]]++;
            nr_of_cells++;
        }
    }

    fprintf(file1, "%g  %d  ", Time, nr_of_cells);
    for(int i=3; i<=12; i++) fprintf(file1, "%d  ", polygonality[i]);
    fprintf(file1, "\n");
    
    delete []polygonality;
    fclose(file1);
}
//****************************************************************************
void aggregateDiameter(){
    
    FILE *file1; file1 = fopen("./output/aggregateDiameter.txt", "at");
    
    int v1, v2, nrEdges=0;
    double avgD=0, avgX, avgY;
    for(int i=1; i<=Ne; i++){
        if(e[i][0]!=0 && e_cell2[i]==0){
            nrEdges++;
            v1=e[i][1]; v2=e[i][2];
            avgX=(v[v1][1]+v[v2][1])/2;
            avgY=(v[v1][2]+v[v2][2])/2;
            avgD+=sqrt(pow(avgX-25,2)+pow(avgY-25,2));
        }
    }
    
    fprintf(file1, "%g  %g\n", Time, avgD/(1.*nrEdges));
    
    fclose(file1);
}
//****************************************************************************
void run(){
    set_initial_fromFile("./initial/hexagonFlat1.vt3d");
    h=0.001;
    int nr_of_cells=Nc, intTime=0;
    double countTime=1000;

    double hess_guess = 1.0;
    bool do_implicit = false;
    bool yolk_present = false;
    while(Time<tmax){
        if(countTime>2){
            countPolygons(); aggregateDiameter();
            intTime++;
            write_vtk(intTime, 0., 0., 0.);
            countTime=0;
        }
        if (eqOfMotion(hess_guess, do_implicit, yolk_present) != 0) {
            printf("Exiting time loop\n");
            break;
        }
        project_on_constraints(30.);
        spont_T1(0.2);
        
        nr_of_cells+=divide_cell(2*h/(1.*nr_of_cells));
        //nr_of_cells+=division_extrusion(0.000013);
        
        countTime+=h;
        if(nr_of_cells>=800) break;
    }
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
