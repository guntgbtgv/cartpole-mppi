#include <eigen3/Eigen/Dense>
#include "custommath.h"
#include <random>
class CController
{

public:

    CController();
    virtual ~CController(){};

    void read(double time, double* q, double* qdot);
    void control_mujoco();
    void write(double* torque);


private:

    VectorXd state_transition(Vector4d state, MatrixXd control, double delta_t);
    double running_cost(VectorXd state, VectorXd control, VectorXd controlvar);
    double terminal_cost(VectorXd state);

    double _t, _dt, _pre_t, _init_t;
    VectorXd _q, _qdot;
    VectorXd _torque;
    bool _bool_init;
    int _n_dof, _a_dof;

    // MPPI related 
    int _K, _N;
    MatrixXd _u, _delta_u;
    double _del_t;
    VectorXd _f;
    MatrixXd _G;
    double _nu;
    double _lambda;
    Matrix2d _R;
    MatrixXd _x;
    VectorXd _S;

    // random sampling from normal distribution
    std::random_device _rd;
    std::mt19937 _gen;
    std::normal_distribution<double> _dist;

};

void dataLogging(double time, double control);
