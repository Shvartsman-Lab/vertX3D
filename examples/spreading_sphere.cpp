//****************GLOBAL PARAMETERS************************
double bet  =1.;    //basal tension
double alph =1.;    //apical tension
double kV   =100;   //reciprocal isothermal compressibility
double V0   =1.;    //preferred volume

//*******************INCLUDE*************************
#include "functions.h"
//***************************************************

void project_on_constraints(double r0){
    //SPHERE
    double S=25., imen, ulomek;
    for(size_t i=1; i<=Nv; i++){
        if(v[i][0]>0.5 && v_type[i]==1){
            imen=(S-v[i][1])*(S-v[i][1])+(S-v[i][2])*(S-v[i][2])+(S-v[i][3])*(S-v[i][3]);
            ulomek=r0*r0/imen;
            v[i][1]+=0.5*(v[i][1]-S)*(ulomek-1.);
            v[i][2]+=0.5*(v[i][2]-S)*(ulomek-1.);
            v[i][3]+=0.5*(v[i][3]-S)*(ulomek-1.);
        }
    }

}
//****************************************************************************
void spont_T1(double _threshold_length){
    for(size_t i = 1; i<=Ne; i++){
        
        if( e[i][0]!=0 ){
            
            double eLen=edge_length(i);
            
            if( eLen<_threshold_length && eLen<e_length[i]){
                std::clog << "Spontaneous ";
                T1(i,5*_threshold_length/4.);
            }
            else e_length[i]=eLen;
        }
    }
}
//****************************************************************************
int division_extrusion(double k){
    int cell_balance=0;
    for(size_t i=1; i<=Nc; i++) if(basal_edges[i][1]!=0) if ( rnd()<k ) cell_balance+=cell_division(i);
    return cell_balance;
}
//****************************************************************************
void run(){
    set_initial_fromFile("./initial/hexagonSphere1.vt3d");
    h=0.001;
    int nr_of_cells=Nc, intTime=0;
    double countTime=1000;

    double hess_guess = 1.0;
    int do_implicit = 0;
    while(Time<tmax){
        if(countTime>2){
            intTime++;
            // outStructure(intTime);
            write_vtk(intTime);
            countTime=0;
        }
        if (eqOfMotion(hess_guess, do_implicit) != 0) {
            printf("Exiting time loop\n");
            break;
        }
        project_on_constraints(5.);
        spont_T1(0.2);
        nr_of_cells+=division_extrusion(2*h/(1.*nr_of_cells));
        countTime+=h;
        if(nr_of_cells>=800) break;
    }
}
//****************************************************************************

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
