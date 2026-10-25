// Independent Section 7/8 check: evaluate the E4 words with GiNaC's G
// functions, not by integrating the complete logarithmic insertion.
// Build: c++ -std=c++17 -O2 checks/ginac_empl_validation.cpp -o checks/ginac_empl_validation $(pkg-config --cflags --libs ginac)
#include <ginac/ginac.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace GiNaC;
using C = std::complex<double>;
using V = std::array<C, 3>;
using M = std::array<std::array<double, 3>, 3>;

struct Rule {
    std::vector<double> x, w;
    explicit Rule(int n) {
        x.resize(n); w.resize(n);
        for (int k=0; k<(n+1)/2; ++k) {
            double t=std::cos(M_PI*(k+0.75)/(n+0.5));
            double p=0, dp=0;
            for (int iter=0; iter<30; ++iter) {
                double a=1, b=0;
                for (int j=1; j<=n; ++j) {
                    double q=((2*j-1)*t*a-(j-1)*b)/j;
                    b=a; a=q;
                }
                p=a; dp=n*(t*a-b)/(t*t-1);
                double dt=p/dp; t-=dt;
                if (std::abs(dt)<3e-16) break;
            }
            x[k]=(1-t)/2; x[n-1-k]=(1+t)/2;
            double wt=1/((1-t*t)*dp*dp);
            w[k]=wt; w[n-1-k]=wt;
        }
    }
};

static C zcomplex(const ex& e) {
    ex v=e.evalf();
    return {ex_to<numeric>(real_part(v).evalf()).to_double(),
            ex_to<numeric>(imag_part(v).evalf()).to_double()};
}

struct Words {
    V e0{}, e1{}, e2{};
    C pair12{}, pair21{}, shuffle{};
};

static Words eval_words(double tau, int n, bool depth_two=false) {
    Rule q(n);
    const double h=std::sqrt(3.0)/2;
    const double uh=std::sqrt(1-h);
    const double eta=1e-12;
    const ex r=std::sqrt(tau);
    const std::array<ex,6> letters={
        numeric(-1), numeric(1), -I*r, I*r,
        -numeric(h)+I*numeric(eta), numeric(h)-I*numeric(eta)};
    const std::array<int,6> signs={1,1,-1,-1,-1,-1};
    const C c4={std::sqrt(tau)/2,0.5};
    Words out;
    // z=1-u^2 removes the square-root endpoint. Splitting at z=h
    // and a sine-squared map smooths the logarithmic branch point.
    for (auto [lo,hi] : {std::pair<double,double>{0,uh},{uh,1}}) {
        for (int k=0; k<n; ++k) {
            const double t=q.x[k];
            const double theta=M_PI*t/2;
            const double u=lo+(hi-lo)*std::pow(std::sin(theta),2);
            const double z=1-u*u;
            const double jac=(hi-lo)*M_PI/2*std::sin(M_PI*t);
            const C psi=c4*2.0/
                (std::sqrt(2-u*u)*std::sqrt(z*z+tau));
            const double weight=q.w[k]*jac;
            C sum=0;
            std::array<C,6> g{};
            for (int j=0; j<6; ++j) {
                // GiNaC's exact one-letter identity G(c;z)=log(1-z/c)
                // avoids a numerical G-series singularity near z=1,h.
                const std::array<C,6> c={
                    C(-1,0), C(1,0), C(0,-std::sqrt(tau)),
                    C(0,std::sqrt(tau)), C(-h,eta), C(h,-eta)};
                g[j]=(j==1)?C(2*std::log(u),0):
                    std::log(C(1,0)-z/c[j]);
                sum+=double(signs[j])*g[j];
            }
            double zp=1;
            for (int i=0; i<3; ++i, zp*=z*z) {
                out.e0[i]+=weight*psi*zp;
                out.e1[i]+=weight*psi*zp*sum;
                out.e2[i]+=weight*psi*zp*sum*sum/2.0;
            }
            if (depth_two) {
                // A genuinely length-three E4 word, in both orders.
                const C g12=zcomplex(G(lst{letters[0],letters[3]},numeric(z)));
                const C g21=zcomplex(G(lst{letters[3],letters[0]},numeric(z)));
                out.pair12+=weight*psi*g12;
                out.pair21+=weight*psi*g21;
                out.shuffle+=weight*psi*g[0]*g[3];
            }
        }
    }
    return out;
}

static std::array<V,3> coefficients(double tau, int n, bool depth_two=false,
                                     Words* words=nullptr) {
    Words e=eval_words(tau,n,depth_two);
    if (words) *words=e;
    const double gamma=zcomplex(Euler).real();
    const C L0={std::log(4.0/15)-2*gamma,-2*M_PI};
    const C norm=2.0/(M_PI*C(std::sqrt(tau)/2,0.5));
    std::array<V,3> v{};
    for (int i=0; i<3; ++i) {
        v[0][i]=norm*e.e0[i];
        v[1][i]=norm*(L0*e.e0[i]+e.e1[i]);
        v[2][i]=norm*((L0*L0/2.0-M_PI*M_PI/3.0)*e.e0[i]
                           +L0*e.e1[i]+e.e2[i]);
    }
    return v;
}

static V matvec(const M& a, const V& v) {
    V b{};
    for (int i=0;i<3;++i)
        for (int j=0;j<3;++j) b[i]+=a[i][j]*v[j];
    return b;
}

static std::pair<M,M> matrices(double t) {
    double g=(1+t)*(3+4*t);
    M a0{{{{-2*t-3.5,2-5.5/t,6/t}},
          {{-1.5,5.5-2*t,-6}},
          {{1.5*t,-9*t-1.5,6*t}}}};
    M a1{{{{0,8/t+4,-4/t}},
          {{0,3/t-1,4}},
          {{0,t-3,3/t+7}}}};
    for (auto* a : {&a0,&a1})
        for (auto& row : *a)
            for (auto& x : row) x/=g;
    return {a0,a1};
}

static void show(const char* name, const C& z) {
    std::cout << name << " = " << z.real() << " " << (z.imag()<0?"-":"+") 
              << " " << std::abs(z.imag()) << " i\n";
}

int main(int argc,char** argv) {
    try {
        Digits=35;
        int n=argc>1?std::stoi(argv[1]):512;
        if (n<8) throw std::runtime_error("quadrature order must be >=8");
        std::cout << std::setprecision(15);
        Words w;
        auto v=coefficients(1,n,true,&w);
        const C notebook[3][3]={
            {{0.8346268416740711,0},{0.38137988175090404,0},
             {0.2782089472246885,0}},
            {{-2.5364577649106788,-4.486630175535494},
             {-1.4147603054479927,-1.7054712099760727},
             {-1.1663331783529634,-1.113662083932178}},
            {{-10.469483771185637,12.408307312980947},
             {-1.78215425289105,5.539008171245812},
             {-0.00972513529943449,4.135943779506695}}};
        std::cout << "tau=mu=s=1; GiNaC E4 words, quadrature order "
                  << n << " per subinterval\n";
        double max_notebook_error=0;
        for (int k=0;k<3;++k)
            for (int i=0;i<3;++i) {
                std::cout << "I" << (2*i+1) << "[" << k << "] ";
                show("GiNaC",v[k][i]);
                max_notebook_error=std::max(max_notebook_error,
                                            std::abs(v[k][i]-notebook[k][i]));
                std::cout << "  |difference from notebook| = "
                          << std::abs(v[k][i]-notebook[k][i]) << "\n";
            }
        show("E4(0,1,1;0,-1,+i;1)",w.pair12);
        show("E4(0,1,1;0,+i,-1;1)",w.pair21);
        std::cout << "shuffle residual = "
                  << std::abs(w.pair12+w.pair21-w.shuffle) << "\n";

        // Numerically evaluate all 36 ordered depth-two MPLs occurring under
        // the E4 outer kernel. This checks the complete epsilon^2 alphabet,
        // including its analytically continued real branch-point letter.
        double all36res=0;
        {
            const double h=std::sqrt(3.0)/2, eta=1e-8, z=0.95;
            const std::array<ex,6> letters={
                numeric(-1),numeric(1),-I,I,
                -numeric(h)+I*numeric(eta),numeric(h)-I*numeric(eta)};
            const std::array<C,6> c={
                C(-1,0),C(1,0),C(0,-1),C(0,1),
                C(-h,eta),C(h,-eta)};
            const std::array<int,6> sign={1,1,-1,-1,-1,-1};
            C s=0, nested=0;
            for (int j=0;j<6;++j) {
                s+=double(sign[j])*std::log(C(1,0)-z/c[j]);
                for (int k=0;k<6;++k)
                    nested+=double(sign[j]*sign[k])*
                        zcomplex(G(lst{letters[j],letters[k]},numeric(z)));
            }
            all36res=std::abs(nested-s*s/2.0);
            std::cout << "all-36-G2 shuffle residual at z=0.95 = "
                      << all36res << "\n";
        }

        // Explicit native GiNaC iterated_integral cross-check at an interior
        // point, where its Taylor series has rapid convergence.
        symbol zz("zz");
        const ex c4=(1+I)/2;
        const ex psi0=user_defined_kernel(c4/sqrt((1-zz*zz)*(1+zz*zz)),zz);
        const C native=zcomplex(iterated_integral(
            lst{psi0,multiple_polylog_kernel(-1)},numeric(1,2),80));
        Rule q(64); C direct=0;
        for (int j=0;j<64;++j) {
            double z=q.x[j]/2;
            direct+=q.w[j]*0.5*C(0.5,0.5)
                    /std::sqrt((1-z*z)*(1+z*z))
                    *zcomplex(G(lst{-1},numeric(z)));
        }
        show("native GiNaC E4(0,1;0,-1;1/2)",native);
        std::cout << "native-vs-quadrature residual = "
                  << std::abs(native-direct) << "\n";

        double max_de_error=0;
        for (double center : {1.0,2.0}) {
            const double dt=0.005;
            auto vc=(center==1.0)?v:coefficients(center,n);
            auto vp=coefficients(center+dt,n);
            auto vm=coefficients(center-dt,n);
            auto vpp=coefficients(center+2*dt,n);
            auto vmm=coefficients(center-2*dt,n);
            auto [a0,a1]=matrices(center);
            for (int k=0;k<3;++k) {
                V rhs=matvec(a0,vc[k]);
                if (k>0) {
                    V b=matvec(a1,vc[k-1]);
                    for (int i=0;i<3;++i) rhs[i]+=b[i];
                }
                double maxres=0, maxrhs=0;
                for (int i=0;i<3;++i) {
                    C lhs=(-vpp[k][i]+8.0*vp[k][i]
                           -8.0*vm[k][i]+vmm[k][i])/(12*dt);
                    maxres=std::max(maxres,std::abs(lhs-rhs[i]));
                    maxrhs=std::max(maxrhs,std::abs(rhs[i]));
                }
                std::cout << "Section 8 DE residual, tau=" << center
                          << ", epsilon order " << k
                          << ": max absolute = " << maxres
                          << ", relative = " << maxres/std::max(1.0,maxrhs)
                          << "\n";
                max_de_error=std::max(max_de_error,maxres);
            }
        }
        bool pass=max_notebook_error<1e-7 && max_de_error<1e-7
            && std::abs(w.pair12+w.pair21-w.shuffle)<1e-10
            && all36res<1e-10 && std::abs(native-direct)<1e-10;
        std::cout << "VALIDATION: " << (pass?"PASS":"FAIL") << "\n";
        return pass?0:2;
    } catch (const std::exception& e) {
        std::cerr << "GiNaC validation failed: " << e.what() << "\n";
        return 1;
    }
}
