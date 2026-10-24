#include "chem/kumo/local_joint_diagonalization.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
int main() {
    try {
        Eigen::Matrix4d r=Eigen::Matrix4d::Identity();
        const double c=std::cos(0.42),s=std::sin(0.42);
        r(0,0)=c;r(1,1)=c;r(0,1)=s;r(1,0)=-s;
        std::vector<Eigen::MatrixXd> gs;
        for(int k=0;k<3;++k) {
            Eigen::Vector4d d;
            d << 1.0+k, 3.0-k*0.3, 5.0+k, 8.0-k;
            gs.push_back((r*d.asDiagonal()*r.transpose()).eval());
        }
        const auto result=chem::kumo::local_joint_diagonalize(gs);
        if(result.final_offdiag>1e-18 || result.final_offdiag>result.initial_offdiag+1e-12 ||
           (result.rotation.transpose()*result.rotation-Eigen::Matrix4d::Identity()).norm()>1e-10)
            throw std::runtime_error("AJD synthetic regression failed");
        std::cout<<"initial="<<result.initial_offdiag<<" final="<<result.final_offdiag
                 <<" sweeps="<<result.sweeps<<" status=OK\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
