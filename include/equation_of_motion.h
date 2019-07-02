#pragma once

#include "functions.h"
#include "LBFGS_helpers.h"

//****************************************************************************
//****************************************************************************
//***********************EQUATION OF MOTION***********************************
//****************************************************************************
//****************************************************************************
void forces(bool yolk_present, double implicit_eps){
    // w = T A
    wA = 0;
    A_tot = 0;
    size_t _nf=Nf;
    #pragma omp parallel for schedule(guided) shared(v, f, v_F) reduction(+:A_tot, wA)
    for(size_t i = 1; i <= _nf; i++){
        if(f[i][0]!=0){
            double _aTot=f_SurfaceTension_force(i, f, implicit_eps);
            wA    += f_T[i]*_aTot;
            A_tot += _aTot;
        }
    }
    
    // w = k_V (V-V_0)^2
    wV=0;
    size_t _nc=Nc;
    #pragma omp parallel for schedule(guided) shared(v_F, v, f) reduction(+:wV)
    for(size_t i = 1; i <= _nc; i++){
        if(basal_edges[i][1]!=0){
            wV += c_VolCompressibility_force(i, f);
        }
    }

    //YOLK
    wY = 0;
    if (yolk_present) {
        wY = YolkCompressibility_force(f);
    }
}
//****************************************************************************
void fix_central_vertices(){
    
    for(size_t i = 1; i<=Nc; i++){
        if(basal_edges[i][1]!=0) {
            //BASAL
            basal_center(i,0);
            //APICAL
            apical_center(i,0);
        }
    }
    
    //LATERAL
    for(size_t i = 1; i<=Ne; i++){
        if(e_v_pass[i]!=0) {
            lateral_center(i,0);
        }
    }
    
}
//****************************************************************************
// void gradient(const Eigen::VectorXd &x, const Eigen::VectorXd &x_prev, Eigen::VectorXd &grad) {
void gradient(LargeVec &x, LargeVec &x_prev, LargeVec &grad, bool yolk_present, double implicit_eps) {
    for (size_t i = 1; i <= Nv; i++) {
        if(v[i][0]>0.5){
            for (size_t j = 1; j <= 3; j++) {
                v[i][j] = x[(i-1)*3 + j-1];
            }
        }
    }
    fix_central_vertices();

    forces(yolk_present, implicit_eps);

    // Eigen::VectorXd force_vec(x.size());
    LargeVec force_vec(x.size());
    for (size_t i = 1; i <= Nv; i++) {
        if(v[i][0]>0.5){
            for (size_t j = 1; j <= 3; j++) {
                // force_vec((i-1)*3 + j-1) = v_F[i][j];
                force_vec[(i-1)*3 + j-1] = v_F[i][j];
            }
        }
    }

    grad = x - x_prev - h * force_vec;

    //RESETS FORCE
    for (size_t i = 1; i <= Nv; i++) {
        if(v[i][0]>0.5){
            for (size_t j = 1; j <= 3; j++) {
                v_F[i][j] = 0;
            }
        }
    }
}
//****************************************************************************
int L_BFGS(double &hess_guess, bool yolk_present, int print, double implicit_eps) {
    int m = 10;
    double tolDiffX = 1e-6;
    double tolDiffGrad = 1e-4;

    int nDof = 3 * Nv;

    LargeVec alpha = LargeVec(m);
    alpha.setZero();
	LargeVec rho = LargeVec(m);
    rho.setZero();
    LargeVec grad(nDof), q(nDof), grad_old(nDof), x_old(nDof), x0(nDof), x_prev(nDof);

    //Move v into x;
    for (size_t i = 1; i <= Nv; i++) {
        if(v[i][0]>0.5){
            for (size_t j = 1; j <= 3; j++) {
                x_prev[(i-1)*3 + j-1] = v[i][j];
            }
        }
    }
    x0 = x_prev;

    gradient(x0, x_prev, grad, yolk_present, implicit_eps);

    VecContainer s(m), y(m);
    for (int i = 0; i < m; i++) {
        s[i].Init(nDof);
        y[i].Init(nDof);
    }

    double gamma_k = hess_guess, gradNorm = 0.0;
	double alpha_init = std::min(1.0, 1.0 / grad.norm());


	int globalIter = 0;

    int k_max = 1000;
    int k = 0;
    bool useGrad = false;
    while (k < k_max) {
        x_old = x0;
        grad_old = grad;
        q = grad;
        globalIter++;

        int iter = std::min(m, k);
		for (int i = iter - 1; i >= 0; --i) {
            rho[i] = 1.0 / s[i].dot(y[i]);
            alpha[i] = rho[i] * s[i].dot(q);
            q -= (alpha[i] * y[i]);
		}

        q = q * gamma_k;

		for (int i = 0; i < iter; i++) {
			double beta = rho[i] * q.dot(y[i]);
			q = q + (alpha[i] - beta) * s[i];
		}

		double dir = q.dot(grad);
        useGrad = false;
		if (dir < 1e-2 * grad.norm() * q.norm()){
            useGrad = true;
			q = grad;
			k_max -= k;
			k = 0;
			alpha_init = std::min(1.0, 1.0 / grad.norm());
		}

        double rate = alpha_init;
		x0 = x0 - rate * q;

        double diffNorm = (x_old - x0).norm() / rate;
		if (diffNorm < tolDiffX){
            if (print > 0) {
                if (useGrad) {
                    printf("t = %f  GD_Iteration: %d    ResNorm = DNC    DiffNorm = %e    rate = %e\n", Time, globalIter, diffNorm, rate);
                } else {
                    printf("t = %f    Iteration: %d    ResNorm = DNC     DiffNorm = %e    rate = %e\n", Time, globalIter, diffNorm, rate);
                }
            }
			break;
		}

		gradient(x0, x_prev, grad, yolk_present, implicit_eps);

		gradNorm = grad.norm();

        if (print > 0) {
            if (useGrad) {
                printf("t = %f  GD_Iteration: %d    ResNorm = %e     DiffNorm = %e    rate = %e\n", Time, globalIter, gradNorm, diffNorm, rate);
            } else {
                printf("t = %f    Iteration: %d    ResNorm = %e     DiffNorm = %e    rate = %e\n", Time, globalIter, gradNorm, diffNorm, rate);
            }
        }

		if (gradNorm < tolDiffGrad){
			hess_guess = gamma_k;
            break;
		}

		if (k < m) {
			s[k] = x0 - x_old;
			y[k] = grad - grad_old;

            gamma_k = s[k].dot(y[k]) / y[k].normSquared();
		} else {
            s.rotate(x0 - x_old);
            y.rotate(grad - grad_old);

            gamma_k = s[m-1].dot(y[m-1]) / y[m-1].normSquared();
		}

		alpha_init = 1.0;
        k++;
    }

    if (k == k_max) {
        return -1;
    }

    for (size_t i = 1; i <= Nv; i++) {
        if(v[i][0]>0.5){
            for (size_t j = 1; j <= 3; j++) {
                v[i][j] = x0[(i-1)*3 + j-1];
            }
        }
    }
    return 0;
}
//****************************************************************************
int eqOfMotion(double &hess_guess, bool do_implicit, bool yolk_present) {
    if (do_implicit)  {
        if (L_BFGS(hess_guess, yolk_present, 1, 1e-3) != 0) {
            printf("Failed to converge\n");
            return -1;
        }
    } else {
        //CALCULATES FORCES
        forces(yolk_present, 0.0);

        //EQUATION OF MOTION
        #pragma omp parallel for schedule(guided) shared(v_F, v)
        for(size_t i = 1; i <= Nv; i++){
            if(v[i][0]>0.5){
                
                double pm_new=0.0;
                
                #pragma omp simd reduction(+:pm_new)
                //MOVES VERTICES
                for(size_t j = 1; j<=3; j++){
                    v[i][j] += h*v_F[i][j];
                    pm_new += pow(v_F[i][j], 2);
                }
                pm_new = h*sqrt(pm_new);
                max_move = pm_new > max_move ? pm_new : max_move;
            
                //PERIODIC BOUNDARY CONDITIONS
                torus_vertex(i);
                
                //RESETS FORCE
                #pragma omp simd
                for(size_t j = 1; j<=3; j++) v_F[i][j] = 0;
            }
        }
    }

    //FIXES PASSIVE VERTICES
    fix_central_vertices();
    
    //t=t+dt
    Time+=h;
    
    printf("%g\tenergyA: %.10f\tenergyV: %.10f\tenergyY: %.10f\tenergy: %.10f\tNv = %d\tNc = %d\tNe = %d\n", Time, wA, wV, wY, wA+wV+wY, Nv, Nc, Ne);

    return 0;
}
//****************************************************************************