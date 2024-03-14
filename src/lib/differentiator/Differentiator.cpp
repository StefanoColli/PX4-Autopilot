#include "Differentiator.hpp"

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <boost/math/special_functions/factorials.hpp>

Differentiator::Differentiator(double r, double Ts) :
    Differentiator(-1.0, r, 1, 0.0, Ts)
{
}

Differentiator::Differentiator(double d, double r, unsigned int m, double mu, double Ts)
{
    // Check parameter values
    if ((d>0.0) || (d<-1.0)) {
        PX4_ERR("d must be a value between -1 and 0.");
    }
    if (r<0.0) {
        PX4_ERR("r must be a positive value.");
    }
    if ((m!=0) && (m!=1) && (m!=2)) {
        PX4_ERR("m must be out of the set {0,1,2}.");
    }
    if (mu<0.0) {
        PX4_ERR("mu must be a positive value.");
    }
    if (Ts<=0.0) {
        PX4_ERR("Ts must be a value greater than 0.");
    }
    if ((m==2) && !(fabs(1 + d) < DBL_EPSILON)) {
        PX4_ERR("Uniform RED (m=2) is only allowed with RED (d=-1).");
    }

    _d  = d;
    _r  = r;
    _mu = mu;
    _Ts = Ts;
    _m  = m;


    // Initialize internal variables
    _zp.setZero();

    _A.setZero();
    //set first super-diagonal to 1.0
    for (unsigned int i = 0; i < _sys_n-1; i++)
    {
        _A(i,i+1) = 1.0;
    }

    _b.setZero();
    if (_nf>0)
    {
        _b(_nf-1) = -1.0;
    }

    _P.setZero();
    _P.slice<_sys_n,_sys_n>(0,0) = _A;
    _P.slice<_sys_n,1>(0,_sys_n) = _b;
    _P = _P*_Ts;

    _P = matrix::exp(_P);

    _Phi = _P.slice<_sys_n,_sys_n>(0,0);
    _bD = _P.slice<_sys_n,1>(0,_sys_n);

    matrix::SquareMatrix<double,_sys_n> fact_diag = matrix::zeros<double, _sys_n, _sys_n>();
    for (unsigned int k=0; k<_sys_n; k++) {
        fact_diag(k,k) = 1/boost::math::factorial<double>(k);
    }

    matrix::Vector<double, _sys_n> ones_sys_n = matrix::ones<double, _sys_n, 1>();
    matrix::Matrix<double, _sys_n, 1> U_vect;
    // lin spaced initialization
    for (unsigned int i = 0; i < _sys_n; i++)
    {
        U_vect(i,0) = i;
    }

    matrix::SquareMatrix<double, _sys_n> U = U_vect * ones_sys_n.transpose();
    matrix::SquareMatrix<double, _sys_n> Ut = U.transpose();

    _U = U.epow(Ut);    //element wise power
    _U = _U*fact_diag;

    _V_inv = matrix::zeros<double, _sys_n, _sys_n>();
    for (unsigned int k=0; k<_sys_n; k++)
    {
        _V_inv(k,k) = 1/std::pow(_Ts,k);
    }

    //_So_inv = (_U.transpose().lu().solve(_V_inv.transpose())).transpose(); (not needed since we use ackerman precomputed)

    _x0 = _z0 = 0.0;
}

Differentiator::~Differentiator()
{
    // Do nothing
}

void Differentiator::evaluate(double u)
{
    _x0 = err(u, _zp(0));

    double z_eig;
    if (_m == 2)
    {
        z_eig = disc_eigenvalues_URED(_x0);
    }
    else
    {
        double s = cont_eigenvalues(_x0);
        z_eig = disc_eigenvalues(s);
    }

    matrix::Vector<double,_sys_n> lambda;
    //if (_sys_n>11) {
    //    ackerman(lambda, z_eig);
    //}
    //else {
    ackerman_precomputed(lambda, z_eig);
    //}

    _z0 = _zp(_nf);
    //last n elements of _zp are the requested n derivatives
    for (unsigned int i = 0; i < _n; i++)
    {
        _z(i) = _zp(_sys_n-_n+i);
    }


    step(_zp, u, lambda);
}


// Private methods
double Differentiator::err(double u, double first_state)
{
    if (_nf>0) {
        return -first_state;
    }

    return u-first_state;
}

double Differentiator::cont_eigenvalues(double x0)
{
    return -_r*std::pow(std::fabs(x0),(double)(_d/(1-_d*(_sys_n-1))));
}

double Differentiator::disc_eigenvalues(double s)
{
    switch (_m) {
        case 0:
            if (std::isinf(s)) {
                return 0.0;
            }
            else {
                return 1.0 + _Ts*s;
            }
            break;

        case 1:
            return std::exp(_Ts*s);
            break;

        default:
            return 0.0;
    }

    return 0.0;
}

double Differentiator::disc_eigenvalues_URED(double x0)
{
    double c = std::pow(std::fabs(x0),1.0/(double)_sys_n);
    return c/(_Ts*_r*_mu*std::abs(x0)+c+_Ts*_r);
}

void Differentiator::ackerman_precomputed(matrix::Vector<double,_sys_n>& lambda, double z)
{
    double z0 = z;

    switch (_sys_n-1) {
        case 1:
            lambda(0) =  2.0-z0-z;
            lambda(1) = ((z-1.0)*(z0-1.0))/_Ts;
            break;

        case 2:
            lambda(0) = 3.0-z0-2.0*z;
            lambda(1) = ((z-1.0)*(z+3.0*z0+z*z0-5.0))/(2.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(z0-1.0))/std::pow(_Ts,2.0);
            break;

        case 3:
            lambda(0) = 4.0-z0-3.0*z;
            lambda(1) = ((z-1.0)*(7.0*z+11.0*z0+5.0*z*z0+2.0*std::pow(z,2.0)*z0+std::pow(z,2.0)-26.0))/(6.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(2.0*z0+z*z0-3.0))/std::pow(_Ts,2.0);
            lambda(3) = (std::pow(z-1.0,3.0)*(z0-1.0))/std::pow(_Ts,3.0);
            break;

        case 4:
            lambda(0) =  5.0-z0-4.0*z,
            lambda(1) = ((z-1.0)*(23.0*z+25.0*z0+13.0*z*z0+7.0*std::pow(z,2.0)*z0+3*std::pow(z,3.0)*z0+5*std::pow(z,2.0)+std::pow(z,3.0)-77.0))/(12.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(35.0*z0-2.0*z+26.0*z*z0+11.0*std::pow(z,2.0)*z0+std::pow(z,2.0)-71.0))/(12.0*std::pow(_Ts,2.0));
            lambda(3) = -(std::pow(z-1.0,3.0)*(z-5.0*z0-3.0*z*z0+7.0))/(2*std::pow(_Ts,3.0));
            lambda(4) = -(std::pow(z-1.0,4.0)*(z0-1.0))/std::pow(_Ts,4.0);
            break;

        case 5:
            lambda(0) = 6.0-z0-5.0*z,
            lambda(1) = ((z-1.0)*(163.0*z+137.0*z0+77.0*z*z0+47.0*std::pow(z,2.0)*z0+27.0*std::pow(z,3.0)*z0+12.0*std::pow(z,4.0)*z0+43.0*std::pow(z,2.0)+13.0*std::pow(z,3.0)+3.0*std::pow(z,4.0)-522.0))/(60.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(45.0*z0-7.0*z+40.0*z*z0+25.0*std::pow(z,2.0)*z0+10.0*std::pow(z,3.0)*z0+2.0*std::pow(z,2.0)+std::pow(z,3.0)-116.0))/(12.0*std::pow(_Ts,2.0));
            lambda(3) = -(std::pow(z-1.0,3.0)*(8.0*z-17.0*z0-16.0*z*z0-7.0*std::pow(z,2.0)*z0+std::pow(z,2.0)+31.0))/(4.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(z-3.0*z0-2.0*z*z0+4.0))/std::pow(_Ts,4.0);
            lambda(5) = (std::pow(z-1.0,5.0)*(z0-1.0))/std::pow(_Ts,5.0);
            break;

        case 6:
            lambda(0) = 7.0-z0-6.0*z,
            lambda(1) = ((z-1.0)*(213.0*z+147.0*z0+87.0*z*z0+57.0*std::pow(z,2.0)*z0+37.0*std::pow(z,3.0)*z0+22.0*std::pow(z,4.0)*z0+10.0*std::pow(z,5.0)*z0+63.0*std::pow(z,2.0)+23.0*std::pow(z,3.0)+8.0*std::pow(z,4.0)+2.0*std::pow(z,5.0)-669.0))/(60.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(812.0*z0-232.0*z+802.0*z*z0+597.0*std::pow(z,2.0)*z0+352.0*std::pow(z,3.0)*z0+137.0*std::pow(z,4.0)*z0+33.0*std::pow(z,2.0)+38.0*std::pow(z,3.0)+13.0*std::pow(z,4.0)-2552.0))/(180.0*std::pow(_Ts,2.0));
            lambda(3) = -(std::pow(z-1.0,3.0)*(39.0*z-49.0*z0-57.0*z*z0-39.0*std::pow(z,2.0)*z0-15.0*std::pow(z,3.0)*z0+9.0*std::pow(z,2.0)+std::pow(z,3.0)+111.0))/(8.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(26.0*z-35.0*z0-38.0*z*z0-17.0*std::pow(z,2.0)*z0+5.0*std::pow(z,2.0)+59.0))/(6.0*std::pow(_Ts,4.0));
            lambda(5) = -(std::pow(z-1.0,5.0)*(3.0*z-7.0*z0-5.0*z*z0+9.0))/(2.0*std::pow(_Ts,5.0));
            lambda(6) = -(std::pow(z-1.0,6.0)*(z0-1))/std::pow(_Ts,6.0);
            break;

        case 7:
            lambda(0) = 8.0-z0-7.0*z,
            lambda(1) = ((z-1.0)*(1851.0*z+1089.0*z0+669.0*z*z0+459.0*std::pow(z,2.0)*z0+319.0*std::pow(z,3.0)*z0+214.0*std::pow(z,4.0)*z0+130.0*std::pow(z,5.0)*z0+60.0*std::pow(z,6.0)*z0+591.0*std::pow(z,2.0)+241.0*std::pow(z,3.0)+101.0*std::pow(z,4.0)+38.0*std::pow(z,5.0)+10.0*std::pow(z,6.0)-5772.0))/(420.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(938.0*z0-414.0*z+994.0*z*z0+819.0*std::pow(z,2.0)*z0+574.0*std::pow(z,3.0)*z0+329.0*std::pow(z,4.0)*z0+126.0*std::pow(z,5.0)*z0+16.0*std::pow(z,2.0)+61.0*std::pow(z,3.0)+36.0*std::pow(z,4.0)+11.0*std::pow(z,5.0)-3490.0))/(180.0*std::pow(_Ts,2.0));
            lambda(3) = -(std::pow(z-1.0,3.0)*(1127.0*z-967.0*z0-1277.0*z*z0-1077.0*std::pow(z,2.0)*z0-647.0*std::pow(z,3.0)*z0-232.0*std::pow(z,4.0)*z0+357.0*std::pow(z,2.0)+77.0*std::pow(z,3.0)+7.0*std::pow(z,4.0)+2632.0))/(120.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(68.0*z-56.0*z0-77.0*z*z0-56.0*std::pow(z,2.0)*z0-21.0*std::pow(z,3.0)*z0+23.0*std::pow(z,2.0)+4.0*std::pow(z,3.0)+115.0))/(6.0*std::pow(_Ts,4.0));
            lambda(5) = -(std::pow(z-1.0,5.0)*(43.0*z-46.0*z0-55.0*z*z0-25.0*std::pow(z,2.0)*z0+10.0*std::pow(z,2.0)+73.0))/(6.0*std::pow(_Ts,5.0));
            lambda(6) = (std::pow(z-1.0,6.0)*(2.0*z-4.0*z0-3.0*z*z0+5.0))/std::pow(_Ts,6.0);
            lambda(7) = (std::pow(z-1.0,7.0)*(z0-1.0))/std::pow(_Ts,7.0);
            break;

        case 8:
            lambda(0) = 9.0-z0-8.0*z,
            lambda(1) = ((z-1.0)*(4437.0*z+2283.0*z0+1443.0*z*z0+1023.0*std::pow(z,2.0)*z0+743.0*std::pow(z,3.0)*z0+533.0*std::pow(z,4.0)*z0+365.0*std::pow(z,5.0)*z0+225.0*std::pow(z,6.0)*z0+105.0*std::pow(z,7.0)*z0+1497.0*std::pow(z,2.0)+657.0*std::pow(z,3.0)+307.0*std::pow(z,4.0)+139.0*std::pow(z,5.0)+55.0*std::pow(z,6.0)+15.0*std::pow(z,7.0)-13827.0))/(840.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(29531.0*z0-18254.0*z+32926.0*z*z0+29013.0*std::pow(z,2.0)*z0+22468.0*std::pow(z,3.0)*z0+15293.0*std::pow(z,4.0)*z0+8622.0*std::pow(z,5.0)*z0+3267.0*std::pow(z,6.0)*z0-733.0*std::pow(z,2.0)+2172.0*std::pow(z,3.0)+1787.0*std::pow(z,4.0)+898.0*std::pow(z,5.0)+261.0*std::pow(z,6.0)-127251.0))/(5040.0*std::pow(_Ts,2.0));
            lambda(3) = -(std::pow(z-1.0,3.0)*(3777.0*z-2403.0*z0-3457.0*z*z0-3302.0*std::pow(z,2.0)*z0-2442.0*std::pow(z,3.0)*z0-1367.0*std::pow(z,4.0)*z0-469.0*std::pow(z,5.0)*z0+1462.0*std::pow(z,2.0)+442.0*std::pow(z,3.0)+87.0*std::pow(z,4.0)+5.0*std::pow(z,5.0)+7667.0))/(240.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(5572.0*z-3207.0*z0-5092.0*z*z0-4682.0*std::pow(z,2.0)*z0-2852.0*std::pow(z,3.0)*z0-967.0*std::pow(z,4.0)*z0+2522.0*std::pow(z,2.0)+772.0*std::pow(z,3.0)+127.0*std::pow(z,4.0)+7807.0))/(240.0*std::pow(_Ts,4.0));
            lambda(5) = -(std::pow(z-1.0,5.0)*(122.0*z-81.0*z0-125.0*z*z0-95.0*std::pow(z,2.0)*z0-35.0*std::pow(z,3.0)*z0+50.0*std::pow(z,2.0)+10.0*std::pow(z,3.0)+154.0))/(6.0*std::pow(_Ts,5.0));
            lambda(6) = (std::pow(z-1.0,6.0)*(42.0*z-39.0*z0-50.0*z*z0-23.0*std::pow(z,2.0)*z0+11.0*std::pow(z,2.0)+59.0))/(4.0*std::pow(_Ts,6.0));
            lambda(7) = -(std::pow(z-1.0,7.0)*(5.0*z-9.0*z0-7.0*z*z0+11.0))/(2.0*std::pow(_Ts,7.0));
            lambda(8) = -(std::pow(z-1.0,8.0)*(z0-1.0))/std::pow(_Ts,8.0);
            break;

        case 9:
            lambda(0) = 10.0-z0-9.0*z,
            lambda(1) = ((z-1.0)*(15551.0*z+7129.0*z0+4609.0*z*z0+3349.0*std::pow(z,2.0)*z0+2509.0*std::pow(z,3.0)*z0+1879.0*std::pow(z,4.0)*z0+1375.0*std::pow(z,5.0)*z0+955.0*std::pow(z,6.0)*z0+595.0*std::pow(z,7.0)*z0+280.0*std::pow(z,8.0)*z0+5471.0*std::pow(z,2.0)+2531.0*std::pow(z,3.0)+1271.0*std::pow(z,4.0)+641.0*std::pow(z,5.0)+305.0*std::pow(z,6.0)+125.0*std::pow(z,7.0)+35.0*std::pow(z,8.0)-48610.0))/(2520.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(32575.0*z0-26477.0*z+37754.0*z*z0+34905.0*std::pow(z,2.0)*z0+28864.0*std::pow(z,3.0)*z0+21689.0*std::pow(z,4.0)*z0+14514.0*std::pow(z,5.0)*z0+8095.0*std::pow(z,6.0)*z0+3044.0*std::pow(z,7.0)*z0-2712.0*std::pow(z,2.0)+2321.0*std::pow(z,3.0)+2566.0*std::pow(z,4.0)+1677.0*std::pow(z,5.0)+788.0*std::pow(z,6.0)+223.0*std::pow(z,7.0)-159826.0))/(5040.0*std::pow(_Ts,2.0));
            lambda(3) = (std::pow(z-1.0,3.0)*(180920.0*z0-363543.0*z+276981.0*z*z0+287607.0*std::pow(z,2.0)*z0+240602.0*std::pow(z,3.0)*z0+165702.0*std::pow(z,4.0)*z0+88737.0*std::pow(z,5.0)*z0+29531.0*std::pow(z,6.0)*z0-161922.0*std::pow(z,2.0)-60422.0*std::pow(z,3.0)-17337.0*std::pow(z,4.0)-2931.0*std::pow(z,5.0)+16.0*std::pow(z,6.0)-663941.0))/(15120.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(9853.0*z-4275.0*z0-7488.0*z*z0-7938.0*std::pow(z,2.0)*z0-6108.0*std::pow(z,3.0)*z0-3363.0*std::pow(z,4.0)*z0-1068.0*std::pow(z,5.0)*z0+5368.0*std::pow(z,2.0)+2198.0*std::pow(z,3.0)+638.0*std::pow(z,4.0)+101.0*std::pow(z,5.0)+12082.0))/(240.0*std::pow(_Ts,4.0));
            lambda(5) = -(std::pow(z-1.0,5.0)*(6428.0*z-3013.0*z0-5444.0*z*z0-5334.0*std::pow(z,2.0)*z0-3284.0*std::pow(z,3.0)*z0-1069.0*std::pow(z,4.0)*z0+3534.0*std::pow(z,2.0)+1244.0*std::pow(z,3.0)+229.0*std::pow(z,4.0)+6709.0))/(144.0*std::pow(_Ts,5.0));
            lambda(6) = (std::pow(z-1.0,6.0)*(129.0*z-75.0*z0-126.0*z*z0-99.0*std::pow(z,2.0)*z0-36.0*std::pow(z,3.0)*z0+60.0*std::pow(z,2.0)+13.0*std::pow(z,3.0)+134.0))/(4.0*std::pow(_Ts,6.0));
            lambda(7) = -(std::pow(z-1.0,7.0)*(172.0*z-145.0*z0-196.0*z*z0-91.0*std::pow(z,2.0)*z0+49.0*std::pow(z,2.0)+211.0))/(12*std::pow(_Ts,7.0));
            lambda(8) = (std::pow(z-1.0,8.0)*(3.0*z-5.0*z0-4.0*z*z0+6.0))/std::pow(_Ts,8.0);
            lambda(9) = (std::pow(z-1.0,9.0)*(z0-1.0))/std::pow(_Ts,9.0);
            break;

        case 10:
            lambda(0) = 11.0-z0-10.0*z,
            lambda(1) = ((z-1.0)*(17819.0*z+7381.0*z0+4861.0*z*z0+3601.0*std::pow(z,2.0)*z0+2761.0*std::pow(z,3.0)*z0+2131.0*std::pow(z,4.0)*z0+1627.0*std::pow(z,5.0)*z0+1207.0*std::pow(z,6.0)*z0+847.0*std::pow(z,7.0)*z0+532.0*std::pow(z,8.0)*z0+252.0*std::pow(z,9.0)*z0+6479.0*std::pow(z,2.0)+3119.0*std::pow(z,3.0)+1649.0*std::pow(z,4.0)+893.0*std::pow(z,5.0)+473.0*std::pow(z,6.0)+233.0*std::pow(z,7.0)+98.0*std::pow(z,8.0)+28.0*std::pow(z,9.0)-55991.0))/(2520.0*_Ts);
            lambda(2) = -(std::pow(z-1.0,2.0)*(177133.0*z0-181196.0*z+211686.0*z*z0+202949.0*std::pow(z,2.0)*z0+175852.0*std::pow(z,3.0)*z0+140985.0*std::pow(z,4.0)*z0+104102.0*std::pow(z,5.0)*z0+68899.0*std::pow(z,6.0)*z0+38136.0*std::pow(z,7.0)*z0+14258.0*std::pow(z,8.0)*z0-27739.0*std::pow(z,2.0)+10278.0*std::pow(z,3.0)+16165.0*std::pow(z,4.0)+12728.0*std::pow(z,5.0)+7611.0*std::pow(z,6.0)+3454.0*std::pow(z,7.0)+962.0*std::pow(z,8.0)-976263.0))/(25200.0*std::pow(_Ts,2.0));
            lambda(3) = (std::pow(z-1.0,3.0)*(420475.0*z0-1040321.0*z+675075.0*z*z0+744585.0*std::pow(z,2.0)*z0+676405.0*std::pow(z,3.0)*z0+526605.0*std::pow(z,4.0)*z0+346845.0*std::pow(z,5.0)*z0+180175.0*std::pow(z,6.0)*z0+58635.0*std::pow(z,7.0)*z0-514467.0*std::pow(z,2.0)-222035.0*std::pow(z,3.0)-80075.0*std::pow(z,4.0)-21303.0*std::pow(z,5.0)-2669.0*std::pow(z,6.0)+427.0*std::pow(z,7.0)-1748357.0))/(30240.0*std::pow(_Ts,3.0));
            lambda(4) = (std::pow(z-1.0,4.0)*(994506.0*z-341693.0*z0-643092.0*z*z0-750990.0*std::pow(z,2.0)*z0-665660.0*std::pow(z,3.0)*z0-462765.0*std::pow(z,4.0)*z0-238632.0*std::pow(z,5.0)*z0-72368.0*std::pow(z,6.0)*z0+617430.0*std::pow(z,2.0)+304040.0*std::pow(z,3.0)+118155.0*std::pow(z,4.0)+33126.0*std::pow(z,5.0)+5084.0*std::pow(z,6.0)+1102859.0))/(15120.0*std::pow(_Ts,4.0));
            lambda(5) = -(std::pow(z-1.0,5.0)*(24135.0*z-8591.0*z0-17305.0*z*z0-19830.0*std::pow(z,2.0)*z0-15730.0*std::pow(z,3.0)*z0-8555.0*std::pow(z,4.0)*z0-2565.0*std::pow(z,5.0)*z0+16010.0*std::pow(z,2.0)+7550.0*std::pow(z,3.0)+2445.0*std::pow(z,4.0)+427.0*std::pow(z,5.0)+22009.0))/(288.0*std::pow(_Ts,5.0));
            lambda(6) = (std::pow(z-1.0,6.0)*(18188.0*z-7513.0*z0-14948.0*z*z0-15378.0*std::pow(z,2.0)*z0-9548.0*std::pow(z,3.0)*z0-3013.0*std::pow(z,4.0)*z0+11418.0*std::pow(z,2.0)+4388.0*std::pow(z,3.0)+853.0*std::pow(z,4.0)+15553.0))/(240.0*std::pow(_Ts,6.0));
            lambda(7) = -(std::pow(z-1.0,7.0)*(1139.0*z-605.0*z0-1085.0*z*z0-875.0*std::pow(z,2.0)*z0-315.0*std::pow(z,3.0)*z0+581.0*std::pow(z,2.0)+133.0*std::pow(z,3.0)+1027.0))/(24.0*std::pow(_Ts,7.0));
            lambda(8) = (std::pow(z-1.0,8.0)*(56.0*z-44.0*z0-62.0*z*z0-29.0*std::pow(z,2.0)*z0+17.0*std::pow(z,2.0)+62.0))/(3.0*std::pow(_Ts,8.0));
            lambda(9) = -(std::pow(z-1.0,9.0)*(7.0*z-11.0*z0-9.0*z*z0+13.0))/(2*std::pow(_Ts,9.0));
            lambda(10) = -(std::pow(z-1.0,10.0)*(z0-1.0))/std::pow(_Ts,10.0);
            break;

        default:
            lambda = matrix::zeros<double, _sys_n, 1>();
            break;
    }
}

void Differentiator::step(matrix::Vector<double, _sys_n>& z, double u, const matrix::Vector<double,_sys_n>& lambda)
{
    z = _Phi*z + _bD*u + lambda*_x0;
}
