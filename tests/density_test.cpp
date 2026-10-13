#include <chem/kumo.hpp>
#include <cmath>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    try {
        for(int f=1;f<3;++f) {
            auto mol=chem::kumo::Molecule::from_fchk(argv[f]);
            auto g=chem::kumo::Grid::around(mol.atoms(),2.0,0.7);
            auto t=mol.density(g,chem::kumo::Spin::Total);
            auto a=mol.density(g,chem::kumo::Spin::Alpha);
            auto b=mol.density(g,chem::kumo::Spin::Beta);
            auto s=mol.density(g,chem::kumo::Spin::Difference);
            for(std::size_t i=0;i<t.values().size();++i) {
                if(!std::isfinite(t.values()[i]) || t.values()[i]<-1e-8 ||
                   std::abs(t.values()[i]-a.values()[i]-b.values()[i])>1e-9 ||
                   std::abs(s.values()[i]-a.values()[i]+b.values()[i])>1e-9) return 1;
            }
            std::cout << argv[f] << " points=" << t.values().size() << " basis=" << mol.basis_size() << '\n';
        }
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
