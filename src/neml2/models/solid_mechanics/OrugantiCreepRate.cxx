// Copyright 2024, UChicago Argonne, LLC
// All Rights Reserved
// Software Name: NEML2 -- the New Engineering material Model Library, version 2
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE (MIT)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "neml2/models/solid_mechanics/OrugantiCreepRate.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/SSR4.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sinh.h"
#include "neml2/tensors/functions/cosh.h"
#include "neml2/tensors/functions/sqrt.h"
//#include "neml2/tensors/functions/dev.h"
//#include "neml2/tensors/functions/norm.h"
//#include "neml2/tensors/functions/outer.h"
//#include "neml2/tensors/functions/imap.h"
#include "neml2/tensors/functions/linalg/eigh.h"

namespace neml2
{
register_NEML2_object(OrugantiCreepRate);

OptionSet
OrugantiCreepRate::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() =
      "Creep rate from Oruganti et al (2011), with multiaxial modification: \\f$ \\dot{\\varepsilon_{ij}} = \\frac{3}{2}\\frac{S_{ij}}{\\sigma_{eq}}\\dot{\\varepsilon^{'}_{0}}"
      "e^{-Q_{C}/RT}\\sinh{\\frac{\\sigma_{eq}(1-H^{*}(1-D_{S}))}{\\sigma_{0}(1-D_{P})}} \\f$"
      "where: ";

  options.set_input("stress") = VariableName(STATE, "S");
  options.set("stress").doc() = "Stress";

  options.set_input("equivalent_stress") = VariableName(STATE, "internal","s");
  options.set("equivalent_stress").doc() = "Equivalent Stress";

  options.set_input("temperature") = VariableName(STATE, "T");
  options.set("temperature").doc() = "Temperature in K";

  options.set_input("mx_evolution") = VariableName(STATE, "internal", "wp");
  options.set("mx_evolution").doc() = "MX Evolution";

  options.set_input("subgrain_evolution") = VariableName(STATE, "internal", "ws");
  options.set("subgrain_evolution").doc() = "Subgrain Evolution";

  options.set_input("isotropic_hardening") = VariableName(STATE, "internal", "k");
  options.set("isotropic_hardening").doc() = "Isotropic Hardening";

  options.set_parameter<TensorName<Scalar>>("stress_scaling");
  options.set("stress_scaling").doc() = "Stress Scaling";

  options.set_parameter<TensorName<Scalar>>("rate_scaling");
  options.set("rate_scaling").doc() = "Rate Scaling";

  options.set_parameter<TensorName<Scalar>>("activation_energy");
  options.set("activation_energy").doc() = "Activation Energy";

  options.set_output("creep_strain_rate") = VariableName(STATE, "internal", "Ec_rate");
  options.set("creep_strain_rate").doc() = "Rate of creep strain";

  return options;
}

OrugantiCreepRate::OrugantiCreepRate(const OptionSet & options)
  : Model(options),
    _S(declare_input_variable<SR2>("stress")),
    _seq(declare_input_variable<Scalar>("equivalent_stress")),
    _T(declare_input_variable<Scalar>("temperature")),
    _wp(declare_input_variable<Scalar>("mx_evolution")),
    _ws(declare_input_variable<Scalar>("subgrain_evolution")),
    _k(declare_input_variable<Scalar>("isotropic_hardening")),
    _hstar(declare_parameter<Scalar>("H", "stress_scaling", true)),
    _edotprime(declare_parameter<Scalar>("e", "rate_scaling", true)),
    _qc(declare_parameter<Scalar>("qc", "activation_energy", true)),
    _Ec_dot(declare_output_variable<SR2>("creep_strain_rate"))
{
}

void
OrugantiCreepRate::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  // Get deviatoric stress
  auto S = SR2(_S).dev();

  auto R = 8.31446261815324; // JK-1mol-1
  auto et = Scalar(exp(-_qc / (R * _T)));
  auto hyp = Scalar((_seq*(1.0-(_hstar*(1.0-_ws)))) / (_k*(1.0-_wp)));
  //auto shyp = sinh(hyp);
  auto dS = SR2::ones_like(SR2(_S))-(1.0/3.0)*SR2::identity(_S.options());


  if (out)
  {
    _Ec_dot =  (3.0 / 2.0) * (S / _seq) * _edotprime * et * sinh(hyp);
  }

  if (dout_din)
  {
    //auto I = imap_v<SR2>(_S.options());
    auto I =SR2::identity(_S.options());
    auto dSds = SR2(dS*(3.0/2.0)*(1.0/_seq) * _edotprime * et * sinh((_seq*(1.0-(_hstar*(1.0-_ws))))/(_k*(1.0-_wp))));

    if (_S.is_dependent()) 
      
      _Ec_dot.d(_S) = dSds.outer(I);

    if (_seq.is_dependent())
      _Ec_dot.d(_seq) = -(3.0/2.0)*_edotprime*S*et*((_k*(_wp-1)*sinh(hyp))+((_hstar*_seq*(_ws-1)+_seq)*cosh(hyp)))/(_k*(_wp-1)*pow(_seq,2.0));
     

    if (_T.is_dependent())
      _Ec_dot.d(_T) = (3.0/2.0)*(S/(_seq*R*pow(_T,2.0))) *_qc* _edotprime * et * sinh(hyp);
  
    if (_wp.is_dependent())
      _Ec_dot.d(_wp) = (3.0/2.0)*_edotprime*S*(_hstar*(_ws-1)+1)*et*cosh(hyp) / (_k*pow(1.0-_wp,2.0)) ;
    
    if (_ws.is_dependent())
      _Ec_dot.d(_ws) = (3.0/2.0)*_edotprime*_hstar*S*et*cosh(hyp) /(_k-_k*_wp);
   
    if (_k.is_dependent())
      _Ec_dot.d(_k) = (3.0/2.0)*_edotprime*S *(_hstar*(_ws-1)+1)*et*cosh(hyp) / (pow(_k,2.0)*(_wp-1)) ;
      
    if (const auto * const hstar = nl_param("hstar"))
      _Ec_dot.d(*hstar) = (-3.0/2.0)*_edotprime*S*(_ws-1)*et*cosh(hyp) / (_k*(_wp-1));

    if (const auto * const edotprime = nl_param("edotprime"))
      _Ec_dot.d(*edotprime) = (3.0/2.0)*(S/_seq) * et * sinh(hyp);

    if (const auto * const qc = nl_param("qc"))
      _Ec_dot.d(*qc) = (-3.0/2.0)*(S/(_seq*R*_T)) * _edotprime * et * sinh(hyp) ;     

  }
}
} // namespace neml2
