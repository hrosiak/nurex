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

#ifndef GlauberModelBase_h
#define GlauberModelBase_h
#include "nurex/Nucleus.h"
#include "nurex/NNCrossSection.h"
#include "nurex/evaporation.h"
namespace nurex{

/**
 * nucleon-nucleon combinations
 */
enum class nucleon_t {pp,pn,nn,np};
enum class nucleonmol_t {target_p, target_n, projectile_p, projectile_n};

/**
 * range types
 */
enum class range_t {ZeroRange, FiniteRange};


/**
 * The enumeration of proton-neutron density combinations for target-projectile
 */
enum special_t {NN=0, pp_Dirac=1, pn_Dirac=2, tp_Dirac=4, tn_Dirac=8, pp_Zero=16, pn_Zero=32, tp_Zero=64, tn_Zero=128};

/**
 * The enumeration of builtin Coulomb corrections
 */
enum class coulomb_correction_t {none=0, classic=1, relativistic=2, sommerfeld=3};

/**
 * The enumeration of builtin Coulomb corrections
 */
enum class cc_correction_t {none=0, PRC82=1, evaporation=2 ,test=3, test2=4 };

struct range_parameters{
  range_parameters(){}
  range_parameters(double beta):pp(beta),pn(beta){}
  range_parameters(double _pp, double _pn):pp(_pp),pn(_pn){}
  bool is_zero()const{return (pp==0.0 && pn==0.0);}
  bool is_same()const{return (pp==pn);}
  double pp=0.0;
  double pn=0.0;
};

using removals_type = std::array<double,6>;

struct EvaporationProbabilities{
    removals_type Ptot{0,0,0,0,0,0};
    removals_type Pch{0,0,0,0,0,0};
    removals_type Pg{0,0,0,0,0,0};
    removals_type Pp{0,0,0,0,0,0};
    removals_type Pn{0,0,0,0,0,0};
    removals_type Pd{0,0,0,0,0,0};
    removals_type Pt{0,0,0,0,0,0};
    removals_type Pa{0,0,0,0,0,0};
    removals_type Phe3{0,0,0,0,0,0};
    removals_type Pimf{0,0,0,0,0,0};
};

/////////////// Glauber Model Base Clas  ////////////////////////
/**
 * Glauber Model Type
 *
 * this class accept any GlauberModel class
 *
*/

class GlauberModelType{
  struct concept_t{
    virtual ~concept_t(){}
    virtual double SigmaR(double)=0;
    virtual double SigmaCC(double)=0;
    virtual double SigmaXN(double)=0;
    virtual double Sigma1N(double)=0;
    virtual EvaporationProbabilities n_removals_evaporation(double)=0;
    virtual removals_type SigmaINs(double)=0;
    virtual void SetCoulombCorrection(coulomb_correction_t)=0;
    virtual void SetCCCorrection(cc_correction_t)=0;
    virtual void SetEvaporationParameters(EvaporationParameters)=0;
    virtual EvaporationParameters GetEvaporationParameters()=0;
    virtual void SetExcitationFunctionType(excitation_function_type par)=0;
    virtual void SetRange(double rpp, double rpn=-1)=0;    
    virtual range_parameters GetRange() const =0;
    virtual const Nucleus& Projectile() const =0;
    virtual const Nucleus& Target() const =0;
    virtual Nucleus ProjectileCopy() const =0;
    virtual Nucleus TargetCopy() const =0;
    virtual void* get_model() =0;
    virtual double X(double, double) = 0;
    virtual double Xpp(double, double) = 0;
    virtual double Xpn(double, double) = 0;
    virtual double Xnp(double, double) = 0;
    virtual double Xnn(double, double) = 0;
  };

  template<typename TYPE>
    struct model_t final:public concept_t{
      static_assert(!std::is_pointer_v<TYPE>, "Shouldn't be pointer here");

      model_t() = default;
      model_t(const TYPE& v) : m_data(v) {}
      model_t(TYPE&& v) : m_data(std::move(v)) {}
      void SetCoulombCorrection(coulomb_correction_t t) override{ m_data.SetCoulombCorrection(t);}
      void SetCCCorrection(cc_correction_t t) override{ m_data.SetCCCorrection(t);}
      void SetEvaporationParameters(EvaporationParameters par) override {m_data.evaporation_parameters = par;}
      EvaporationParameters GetEvaporationParameters() override {return m_data.evaporation_parameters;};
      void SetExcitationFunctionType(excitation_function_type par) override {m_data.SetExcitationFunctionType(par);}
      double SigmaR(double E) override {return m_data.SigmaR(E);}
      double SigmaCC(double E) override {return m_data.SigmaCC(E);}
      double SigmaXN(double E) override {return m_data.SigmaXN(E);}
      double Sigma1N(double E) override {return m_data.Sigma1N(E);}
      removals_type SigmaINs(double E) override {return m_data.SigmaINs(E);}
      EvaporationProbabilities n_removals_evaporation(double E) override {return m_data.n_removals_evaporation(E);}
      void SetRange(double rpp, double rpn=-1) override {m_data.SetRange(rpp, rpn);}
      range_parameters GetRange() const override {return m_data._range;};
      const Nucleus& Projectile() const override{return m_data.Projectile();}
      const Nucleus& Target() const override{return m_data.Target();}
      Nucleus ProjectileCopy() const override{return m_data.ProjectileCopy();}
      Nucleus TargetCopy() const override{return m_data.TargetCopy();}
      void* get_model() override {return &m_data;}
      double X(double b, double E)override{return m_data.X(b, E);}
      double Xpp(double b, double E)override{if constexpr(can_do_cc<TYPE>)return m_data.Xpp(b, E);else return -1.0;}
      double Xpn(double b, double E)override{if constexpr(can_do_cc<TYPE>)return m_data.Xpn(b, E);else return -1.0;}
      double Xnp(double b, double E)override{if constexpr(can_do_cc<TYPE>)return m_data.Xnp(b, E);else return -1.0;}
      double Xnn(double b, double E)override{if constexpr(can_do_cc<TYPE>)return m_data.Xnn(b, E);else return -1.0;}
    TYPE m_data;
  };

  private:
        std::unique_ptr<concept_t> object;

public:
    GlauberModelType()=default;
    GlauberModelType(const GlauberModelType&)=delete;
    GlauberModelType(GlauberModelType&&) = default;
    GlauberModelType& operator=(const GlauberModelType&) = delete;
    GlauberModelType& operator=(GlauberModelType&&) = default;

    template<typename T>
    GlauberModelType(T x):object(std::make_unique<model_t<T>>(std::move(x))){}

    void SetCoulombCorrection(coulomb_correction_t t){object->SetCoulombCorrection(t);}
    void SetCCCorrection(cc_correction_t t){object->SetCCCorrection(t);}
    void SetEvaporationParameters(EvaporationParameters par){object->SetEvaporationParameters(par);}
    void SetExcitationFunctionType(excitation_function_type par) {object->SetExcitationFunctionType(par);}
    EvaporationParameters GetEvaporationParameters(){return object->GetEvaporationParameters();}

    double SigmaR(double E){return object->SigmaR(E);}
    double SigmaCC(double E){return object->SigmaCC(E);}
    double SigmaXN(double E){return object->SigmaXN(E);}
    double Sigma1N(double E){return object->Sigma1N(E);}
    removals_type SigmaINs(double E){return object->SigmaINs(E);}
    EvaporationProbabilities n_removals_evaporation(double E) {return object->n_removals_evaporation(E);}
    void SetRange(double rpp, double rpn=-1) {object->SetRange(rpp, rpn);}
    range_parameters GetRange() const {return object->GetRange();}    
    Nucleus Projectile() const {return object->Projectile();}
    Nucleus Target() const {return object->Target();}
    void* get_object(){return  object->get_model();}
    explicit operator bool()const{return (object)?true:false;}
    double X(double b, double E){return object->X(b,E);}
    double Xpp(double b, double E){return object->Xpp(b,E);}
    double Xpn(double b, double E){return object->Xpn(b,E);}
    double Xnp(double b, double E){return object->Xnp(b,E);}
    double Xnn(double b, double E){return object->Xnn(b,E);}
    double Xcc(double b, double E){return object->Xpp(b,E)+object->Xpn(b,E);}
};


} // namespace nurex
#endif
