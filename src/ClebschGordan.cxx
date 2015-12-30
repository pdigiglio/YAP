#include "ClebschGordan.h"

#include "Exceptions.h"

#include <algorithm>
#include <cmath>

namespace yap {

//-------------------------
bool ClebschGordan::nonzeroClebschGordan(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    // check input spin-projection compatibilities
    if (!consistent(two_j1, two_m1) or !consistent(two_j2, two_m2) or !consistent(two_J,  two_M ))
        throw exceptions::InconsistentSpinProjection();
    // and that (j1+j2) and J are consistent
    if (is_odd(two_J + two_j1 + two_j2))
        throw exceptions::AngularMomentumNotConserved();

    // check input spin-projections
    if (two_M != two_m1 + two_m2)
        return false;

    // check whether J lies between |j1 - j2| and (j1 + j2)
    if (two_J < abs((int)two_j1 - (int)two_j2) or two_J > (two_j1 + two_j2))
        return false;

    // when either daughter spin is zero
    if (two_j1 == 0 and two_j2 != two_J)
        return false;
    if (two_j2 == 0 and two_j1 != two_J)
        return false;

    // check for (j1 0 j2 0 | J 0), j1 + j2 + J is odd
    if (two_m1 == 0 and two_m2 == 0 and is_odd((two_j1 + two_j2 + two_J) / 2))
        return false;

    /// \todo find formulae that apply for these three cases
    /// @{

    // (3/2 +-1/2, 3/2 +-1/2 | 2 +-1) == 0
    if (two_j1 == 3 and abs(two_m1) == 1 and two_j2 == two_j1 and two_m2 == two_m1 and two_J == 4)
        return false;

    // (2 +-1, 3/2 -+1/2 | 3/2 +-1/2) == 0
    if (two_j1 == 4 and abs(two_m1) == 2 and two_j2 == 3 and abs(two_m2) == 1 and two_J == 3 and abs(two_M) == 1)
        return false;

    // (2 +-1, 2 +-1 | 3 +-2) == 0
    if (two_j1 == 4 and abs(two_m1) == 2 and two_j2 == two_j1 and two_m2 == two_m1 and two_J == 6)
        return false;

    /// @}

    return true;
}


//-------------------------
double ClebschGordan::coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    if (!nonzeroClebschGordan(two_j1, two_m1, two_j2, two_m2, two_J, two_M))
        return 0;

    // z range dictated by factorials in denominator ( 1/n! = 0 when n < 0)
    unsigned z_min = std::max({0, (int)two_j2 - two_m1 - (int)two_J, (int)two_j1 + two_m2 - (int)two_J}) / 2;
    unsigned z_max = std::min({two_j1 + two_j2 - two_J, two_j1 - two_m1, two_j2 + two_m2}) / 2;

    double z_sum = 0;
    for (unsigned z = z_min; z <= z_max; ++z) {
        // z'th term := (-)^z / z! / (j1+j2-J-z)! / (j1-m1-z)! / (j2+m2-z)! / (J-j2+m1+z)! / (J-j1-m2+z)!
        z_sum += pow_negative_one(z)
                 / std::tgamma(1 + z)
                 / std::tgamma(1 + (two_j1 + two_j2 - two_J) / 2 - z)
                 / std::tgamma(1 + (two_j1 - two_m1) / 2 - z)
                 / std::tgamma(1 + (two_j2 + two_m2) / 2 - z)
                 / std::tgamma(1 + (two_J - two_j2 + two_m1) / 2 + z)
                 / std::tgamma(1 + (two_J - two_j1 - two_m2) / 2 + z);
    }

    // C-G coef = sqrt( (2J+1)! (j1+j2-J)! (j1-j2+J)! (j2-j1+J)! / (J+j1+j2+1)! )
    //          * sqrt( (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (J+M)! (J-M)! )
    //          * z_sum
    return z_sum * sqrt(std::tgamma(1 + two_J + 1)
                        * std::tgamma(1 + (two_j1 + two_j2 - two_J) / 2)
                        * std::tgamma(1 + (two_j1 - two_j2 + two_J) / 2)
                        * std::tgamma(1 + (two_j2 - two_j1 + two_J) / 2)
                        / std::tgamma(1 + (two_j1 + two_j2 + two_J) / 2 + 1)
                        * std::tgamma(1 + (two_j1 + two_m1) / 2)
                        * std::tgamma(1 + (two_j1 - two_m1) / 2)
                        * std::tgamma(1 + (two_j2 + two_m2) / 2)
                        * std::tgamma(1 + (two_j2 - two_m2) / 2)
                        * std::tgamma(1 + (two_J + two_M) / 2)
                        * std::tgamma(1 + (two_J - two_M) / 2));
}

//-------------------------
double ClebschGordan::couple(unsigned two_j1, int two_lambda1, unsigned two_j2, int two_lambda2, unsigned l, unsigned two_s, unsigned two_J)
{
    if (!nonzeroCoupling(two_j1, two_lambda1, two_j2, two_lambda2, l, two_s, two_J))
        return 0;

    // calculate C-G coefficient for l-s coupling
    double l_s_coupling = coefficient(2 * l, 0, two_s, two_lambda1 - two_lambda2, two_J);
    // calculate C-G coefficient for s-s coupling
    double s_s_coupling = coefficient(two_j1, two_lambda1, two_j2, -two_lambda2, two_s);

    // sqrt( (2l + 1) / (2J + 1) ) * l_s_coupling * s_s_coupling
    return sqrt((2 * l + 1) / (two_J + 1)) * l_s_coupling * s_s_coupling;
}

}
