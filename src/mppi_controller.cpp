#include "mppi_controller.h"
#include <iostream>
#include <math.h>
#include <fstream>

#define mp 0.1
#define mc 1.0
#define l 0.5
#define g 9.81


std::ofstream outFile_mppi("mppi.csv");
bool isFirstloop = true;

CController::CController(){

    _n_dof = 2;
	_a_dof = 1;
    _bool_init = true;
	_t = 0.0;
	_init_t = 0.0;
	_pre_t = 0.0;
	_dt = 0.0;

    _torque.setZero(_a_dof);

    _q.setZero(_n_dof);
	_qdot.setZero(_n_dof);




	_K = 100;
	_N = 10;
	_u.setZero(2,_N);
	// _u.block(0,0,1,_N).setOnes();
	_delta_u.setZero(2*_K,_N);
	_f.setZero(4);
	_G.setZero(4,2);
	_del_t = 0.1;
	_x.setZero(4,_N+1);
	_S.setZero(_K);
	_nu = pow(25.0, 2);
	_R << 1.0, 0.0,
		  0.0, 1.0;
	_lambda = 200.0;
	

	_gen.seed(_rd());
	_dist = std::normal_distribution<double>(0.0, sqrt(_nu));


}


void CController::read(double t, double* q, double* qdot)
{	
	_t = t;
	if (_bool_init == true)
	{
		_init_t = _t;
		_bool_init = false;
	}

	_dt = t - _pre_t;
	// cout<<"_dt : "<<_dt<<endl;
	_pre_t = t;

	for (int i = 0; i < _n_dof; i++)
	{
		_q(i) = q[i];
		_qdot(i) = qdot[i];
	}
}

void CController::control_mujoco(){


	for(int k = 0; k < _K; k++){

		_x.block(0,0,4,1) << _q(0), _q(1), _qdot(0), _qdot(1);
		// _S(k) = running_cost(_x.block(0,0,4,1), _u.block(0,0,2,1), _delta_u.block(2*k,0,2,1));
		// _delta_u.block(2*k,0,2,1) << _dist(_gen), 0.0;
		
		for(int i = 0; i < _N; i++){
			
			_delta_u.block(2*k,i,2,1) << _dist(_gen), 0.0;

			_x.block(0,i+1,4,1) = _x.block(0,i,4,1) + state_transition(_x.block(0,i,4,1), _u.block(0,i,2,1) + _delta_u.block(2*k,i,2,1), _del_t);
			_S(k) += running_cost(_x.block(0,i,4,1), _u.block(0,i,2,1), _delta_u.block(2*k,i,2,1));
		}

		_S(k) += terminal_cost(_x.block(0,_N,4,1));

	}
	// std::cout << "_S(50): " << _S(50) << "\n";
	// std::cout << "_x: \n" << _x << "\n"; 


	long double denominator = 0.0;
	for(int k=0;k<_K;k++){

		denominator += exp((-1.0/_lambda)*_S(k));

	}

	std::cout << "denom: " << denominator << "\n";

	for(int i=0; i<_N; i++){

		for(int k=0;k<_K;k++){

			_u(0,i) += ((exp((-1.0/_lambda)*_S(k)))/denominator)*_delta_u(2*k,i);

		}
	}
	// std::cout << "_u: \n" << _u << "\n";


	_torque(0) = _u(0,0);

	for(int i = 0; i < _N-1; i++){
		_u(0,i) = _u(0,i+1);
	}

	_S.setZero();
	_x.setZero();

	dataLogging(_t, _torque(0));

					// Vector4d K;
					// K << -5.0, 500.0, -10.0, 20.0;

					// Vector4d x, x_d;
					// x << _q(0), _q(1), _qdot(0), _qdot(1); 
					// x_d << 0., 3.14159265358979, 0., 0.; 

					// _torque(0) = -K.dot((x - x_d));

}


VectorXd CController::state_transition(Vector4d state, MatrixXd control, double delta_t){

	
	double theta = state(1);
	double xdot = state(2);
	double thetadot = state(3);
	
	_f(0) = xdot;
	_f(1) = thetadot;
	_f(2) = 1/(mc + mp*sin(theta)*sin(theta)) * (mp*sin(theta)*(l*thetadot*thetadot + g*cos(theta)));
	_f(3) = 1/(l*(mc + mp*sin(theta)*sin(theta))) * (-mp*l*thetadot*thetadot*cos(theta)*sin(theta) - (mc+mp)*g*sin(theta));

	_G.block(2,0,2,2) << 1/(mc + mp*sin(theta)*sin(theta)), 1.0,
						 -cos(theta)/(l*(mc + mp*sin(theta)*sin(theta))), 1.0; 

	
	return (_f + _G*control)*delta_t;

	

}

double CController::running_cost(VectorXd state, VectorXd control, VectorXd controlvar){
	double x = state(0); 
	double theta = state(1); 
	double xdot = state(2);
	double thetadot = state(3);

	double q;

	q = 10.0*x*x + 500.0*pow((1+cos(theta)), 2) + thetadot*thetadot + 10.0*xdot*xdot;

	return q + 0.5*(1.0 - 1.0/_nu)*controlvar.transpose()*_R*controlvar + control.transpose()*_R*controlvar + 0.5*control.transpose()*_R*control;
}

double CController::terminal_cost(VectorXd state){
	double x = state(0); 
	double theta = state(1); 
	double xdot = state(2);
	double thetadot = state(3);

	double q;

	q = 10.0*x*x + 500.0*pow((1+cos(theta)), 2) + thetadot*thetadot + 10.0*xdot*xdot;

	return q;
}



void CController::write(double* torque){

	for (int i = 0; i < _a_dof; i++)
	{
		torque[i] = _torque(i);
	}
}

void dataLogging(double time, double control){

	if (isFirstloop) {
	outFile_mppi << "Time";
	outFile_mppi << ",u";
	outFile_mppi << "\n" ;


	isFirstloop = false;
	}

	outFile_mppi << time;
	outFile_mppi << "," << control;
	outFile_mppi << "\n";  


}