#ifndef DIFFERENZIATOR_H_
#define DIFFERENZIATOR_H_

#include <mathlib/mathlib.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Differentiator class
//
// Implements the differentiator described in
//
// B. Andritsch, M. Horn, S. Koch, H. Niederwieser, M. Wetzlinger, M. Reichhartinger
// The Robust Exact Differentiator Toolbox revisited: Filtering and Discretization Features
// 2021 IEEE International Conference on Mechatronics
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Differentiator
{
    private:
        //diff parameters (treated as constexpr since PX4 mathlib doesn't support dynamic matrices)
        static constexpr unsigned int _n = 1;   // differentiation order (we need jerk)
        static constexpr unsigned int _nf = 2;  // filter order

        static constexpr unsigned int _sys_n = _n + _nf + 1;

        double _d, _r, _mu, _Ts;
        unsigned int _m;
        double _x0, _z0;

        matrix::Vector<double, _sys_n> _zp, _b, _bD;
        matrix::Vector<double, _n> _z;
        matrix::SquareMatrix<double, _sys_n> _A,  _Phi, _V_inv, _U, _So_inv;
        matrix::SquareMatrix<double, _sys_n +1> _P;

        double err(double u, double first_state);
        double cont_eigenvalues(double x0);
        double disc_eigenvalues(double s);
        double disc_eigenvalues_URED(double x0);

        void   ackerman_precomputed(matrix::Vector<double,_sys_n>& lambda, double z);
        void   step(matrix::Vector<double, _sys_n>& z, double u, const matrix::Vector<double,_sys_n>& lambda);

    public:
        Differentiator(double r, double Ts);
        Differentiator(double d, double r, unsigned int m, double mu, double Ts);
        ~Differentiator();

        void evaluate(double u);

        void get_x0(double& x0) { x0 = _x0; };
        void get_z0(double& z0) { z0 = _z0; };
        void get_z(matrix::Vector<double,_n>& z) { z = _z; };
};

class Linear_differentiator : public Differentiator
{
    public:
        Linear_differentiator(double r, unsigned int m, double Ts) : Differentiator(0.0, r, m, 0.0, Ts) {};
        ~Linear_differentiator() {};
};


class Robust_exact_differentiator : public Differentiator
{
    public:
        Robust_exact_differentiator(double r, unsigned int m, double Ts) : Differentiator(-1.0, r, m, 0.0, Ts) {};
        ~Robust_exact_differentiator() {};
};


class Uniform_robust_exact_differentiator : public Differentiator
{
    public:
        Uniform_robust_exact_differentiator(double r, double mu, double Ts) : Differentiator(-1.0, r, 2.0, mu, Ts) {};
        ~Uniform_robust_exact_differentiator() {};
};

#endif // DIFFERENZIATOR_H_
