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
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sinh.h"
#include "neml2/tensors/functions/cosh.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/dev.h"
#include "neml2/tensors/functions/norm.h"
#include "neml2/tensors/functions/linalg/eigh.h"

namespace neml2
{
register_NEML2_object(OrugantiCreepRate);

OptionSet
OrugantiCreepRate::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() =
      "Oruganti creep rate. DOC TBC for use with Liu-Murakami Damage Rate."
      "Creep rate given by \\f$ \\dot{\\varepsilon}_{c} = \\frac{3}{2} B \\sigma_{EQ}^{n-1} s_{ij}"
      "\\left[ \\frac{2(n+1)}{\\pi \\sqrt(1+3/n)} \\left) \\frac{\\sigma_{I}}{\\sigma_{EQ}} "
      "\\right)^{2}"
      "\\omega^{3/2} where: ";

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
  auto S = neml2::dev(_S());

  const auto R = 8.31446261815324; // JK-1mol-1
  auto et = exp(-_qc / (R * _T()));

  if (out)
  {
    _Ec_dot = SR2((3.0/2.0)*(S/_seq) * _edotprime * et * sinh((_seq*(1.0-(_hstar*(1.0-_ws))))/(_k*(1.0-_wp))),_S.dynamic_sizes(), _S.intmd_dim() );
  }

  if (dout_din)
  {
    // if (_S.is_dependent()) //Right now not sure how to handle stress dependence
    // //as damage stress derivative is a pain with s1 eigenvalue
    //   _Ec_dot.d(_S) = SR2::identity_map(S.options()) *
    //   Scalar((3.0/2.0)*_B*pow(vm,_n-1)*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)));

    // if (_w.is_dependent())
    //   _Ec_dot.d(_w) =
    //   SR2((9*_B*(_n+1)*S*pow(s1,2.0)*sqrt(_w)*pow(vm,_n-3.0)*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)))/(2*pi*sqrt((_n+3/_n))));

    // if (const auto * const B = nl_param("B"))
    //   _Ec_dot.d(*B) =
    //   SR2((3.0/2.0)*pow(vm,_n-1)*S*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)));

    // if (const auto * const n = nl_param("n"))
    //   _Ec_dot.d(*n) =
    //   SR2(3.0*_B*S*pow(vm,_n-3)*exp(((2*(_n+1))/(pi*sqrt(1+(3/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0))*(pi*pow(_n,2.0)*pow((_n+3)/_n,3.0/2.0)*pow(vm,2.0)*log(vm)
    //   + (2*pow(_n,2.0)+9*_n +3)*pow(s1,2.0)*pow(_w,3.0/2.0) ) );

    // if (const auto * const q = nl_param("q"))
    //    _w_dot.d(*q) = (_A/pow(_q,2.0))*pow(sd,_p)*exp(_q*(_w-1))*(exp(_q)*(_q*_w -1)+_q + 1);

    // if (const auto * const alpha = nl_param("alpha"))
    //    _w_dot.d(*alpha) = (_A/_q)*_p*(exp(_q)-1)*exp(_q*(_w-1))*(vm-sd)*pow(sd,_p-1);

    //   }
  }
}
} // namespace neml2
