
#include "iiwa14.h"

#define DoP 4 // degree of parallelism 

pinocchio::Model model;

// initialize your own size 
// data for the NMPC controller
pinocchio::Data data[DoP] = {pinocchio::Data(model),pinocchio::Data(model),pinocchio::Data(model),pinocchio::Data(model)};

// data for the simulation model
pinocchio::Data sim_data(model);

void iiwa14_init()
{
    std::string filename = "../iiwa_description/iiwa14.urdf";
    pinocchio::urdf::buildModel(filename,model);
    model.gravity.setZero();
    // init for NMPC
    for(int i=0;i<DoP;i++)
    {
        data[i] = pinocchio::Data(model);
    }
    // init for Simulation
    sim_data = pinocchio::Data(model);

}

void qdd_cal(const double *q, const double *qd, double *qdd, const double *tau, const int parIdx)
{
    const ConstMapVector q_Eigen = ConstMapVector(q, model.nv);
    const ConstMapVector qd_Eigen = ConstMapVector(qd,model.nv);
    const ConstMapVector tau_Eigen = ConstMapVector(tau,model.nv);
    MapVector qdd_Eigen = MapVector(qdd, model.nv);
    
    qdd_Eigen = pinocchio::aba(model,data[parIdx-1],q_Eigen,qd_Eigen,tau_Eigen);
}

void sim_qdd_cal(const double *q, const double *qd, double *qdd, const double *tau)
{
    const ConstMapVector q_Eigen = ConstMapVector(q, model.nv);
    const ConstMapVector qd_Eigen = ConstMapVector(qd,model.nv);
    const ConstMapVector tau_Eigen = ConstMapVector(tau,model.nv);
    MapVector qdd_Eigen = MapVector(qdd, model.nv);
    
    qdd_Eigen = pinocchio::aba(model,sim_data,q_Eigen,qd_Eigen,tau_Eigen);
}

void derivatives_cal(const double *q, const double *qd, const double *tau, double *dq, double *dqd, double *dtau, const int parIdx)
{
    const ConstMapVector q_Eigen = ConstMapVector(q, model.nv);
    const ConstMapVector qd_Eigen = ConstMapVector(qd,model.nv);
    const ConstMapVector tau_Eigen = ConstMapVector(tau,model.nv);

    MapMatrix dq_Eigen   = MapMatrix(dq,  model.nv,model.nv);
    MapMatrix dqd_Eigen  = MapMatrix(dqd, model.nv,model.nv);
    MapMatrix dtau_Eigen = MapMatrix(dtau,model.nv,model.nv);
    
    computeABADerivatives(model, data[parIdx-1], q_Eigen, qd_Eigen, tau_Eigen, dq_Eigen, dqd_Eigen, dtau_Eigen);
}
