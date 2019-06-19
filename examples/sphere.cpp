//****************INPUT PARAMETERS************************
double bet  =.5;    //basal tension
double alph =.5;    //apical tension
double kV   =100.;  //reciprocal isothermal compressibility
double V0   =1.;    //preferred volume

// double VY0  =17804.44;
// double VY0 = 3891.4143;

// double kT1=0;

//*******************DECLARATIONS*************************
#include "functions.h"
//********************************************************

//********************************************************
void run(){
    set_initial_fromFile("./initial/out_Final.vt3d");
    h=0.01;
    int intTime=0;
    double countTime=1000;

    VY0 = 19519.0193;
    c_kY =100.0;

    double hess_guess = 1.0;
    int do_implicit = 1;
    bool yolk_present = true;
    while (Time < tmax) {
        if (countTime > 1) {
            intTime++;
            write_vtk(intTime);
            countTime=0;
        }

        if (eqOfMotion(hess_guess, do_implicit, yolk_present) != 0) {
            printf("Exiting time loop\n");
            break;
        }

        countTime+=h;
    }
}

//********************MAIN********************************
int main(int argc, char *argv[]) {
    read_code_arguments(argc,argv);
    allocate();
    run();
    out_Vertissue3D("./output/out_Final.vt3d");
    deallocate();
    return 0;
}