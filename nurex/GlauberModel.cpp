/*
 *  Copyright(C) 2017, Andrej Prochazka, M. Takechi
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "nurex/GlauberModel.h"
#include "nurex/Models.h"
#include "nurex/epax.h"
namespace nurex {

template class GlauberModel<OLA>;
template class GlauberModel<MOL>;
template class GlauberModel<MOL4C>;
template class GlauberModel<OLA, FermiMotion<NNCrossSectionFit>>;
template class GlauberModel<MOL, FermiMotion<NNCrossSectionFit>>;
template class GlauberModel<MOL4C, FermiMotion<NNCrossSectionFit>>;

double coulomb_correction_simple(const Nucleus& projectile, const Nucleus& target, double E, double cs)
{
    double r = 0;
    if(cs<=0.0 || E<=0.0)return 0.0;
    double bmax = sqrt(cs / 10.0 / PI);
    double Bc = 1.44 * projectile.Z() * target.Z() / bmax;
    double Ecm = Ecm_from_T(E,projectile.A(),target.A());
    r = 1.0 - (Bc / Ecm);
    if (r < 0){
        r = 0.0;
        }
    return r;
}

double coulomb_correction_relativistic(const Nucleus& projectile, const Nucleus& target, double E, double cs)
{
    double r = 0;
    if(cs<=0.0 || E<=0.0)return 0.0;
    double bmax = sqrt(cs / 10.0 / PI);
    double Bc = 1.44 * projectile.Z() * target.Z() / bmax;
    double Ecm = Ecm_from_T_relativistic(E,projectile.A(),target.A());
    assert(Ecm>0.0);
    r = 1.0 - (Bc / Ecm);
    if (r < 0){
        r = 0.0;
        }
    return r;
}

double cc_evaporation_cor(Nucleus const &Projectile, const removals_type &nrcs, const EvaporationParameters &par){
  double cor = 0.0;
  double Ex = Emax(Projectile, par);
  
  for(int i=0;i<nrcs.size();i++){
    if(nrcs[i]<1.0)continue; // if cross-section less than 1.0 ignore it, below req. precision
    double P = ( (par.config&evaporation_config_type::disable_neutron) == 0)?
              charge_evaporation_probability_total(Projectile.A()-i-1, Projectile.Z(), Ex, i+1, par)
              :charge_evaporation_probability_simple(Projectile.A()-i-1, Projectile.Z(), Ex, i+1);
    cor += P*nrcs[i];
  }
  return cor;
};

double total_evaporation_cor(Nucleus const &Projectile, double nrcs, const EvaporationParameters &par){
  double cor = 0.0;
  double Ex = Emax(Projectile, par);

  if(nrcs<1.0)return nrcs; // if cross-section less than 0.5 ignore it, below req. precision
  double P = total_evaporation_probability(Projectile.A()-1, Projectile.Z(), Ex, 1, par);
  return P*nrcs;
};

removals_type epax_xn_ratios(const Nucleus&  Projectile, const Nucleus& Target, const double norm){
    removals_type n;
    double sum=0.0;

    for(int i=0;i<n.size();i++){
          n[i] = epax::epax3(Projectile.A(),Projectile.Z(),Target.A(),Target.Z(),Projectile.A()-i-1,Projectile.Z());
          sum+=n[i];
    }

    for(int i=0;i<n.size();i++){
        n[i] *= norm/sum;
    }
    return n;
};


} // namespace nurex
