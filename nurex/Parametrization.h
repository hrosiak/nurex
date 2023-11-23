#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
namespace nurex{
class Nucleus;
double SigmaR_Kox(int Ap, int Zp, double E, int At, int Zt);
double SigmaR_Kox(const Nucleus &projectile, double E, const Nucleus &target);

}
#endif
