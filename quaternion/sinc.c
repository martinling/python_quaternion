#include <math.h>
#include <float.h>

double
sinc(double x)
{
    const double taylor_0_bound = DBL_EPSILON;
    const double taylor_2_bound = sqrt(taylor_0_bound);
    const double taylor_n_bound = sqrt(taylor_2_bound);
    double abs_x = fabs(x);
    if (abs_x >= taylor_n_bound) {
        return sin(x) / x;
    } else {
        double result = 1;
        if (abs_x >= taylor_0_bound)
        {
            double x_squared = x * x;
            result -= x_squared / 6;
            if (abs_x >= taylor_2_bound)
            {
                result += (x_squared * x_squared) / 120;
            }
        }
        return result;
    }
}
