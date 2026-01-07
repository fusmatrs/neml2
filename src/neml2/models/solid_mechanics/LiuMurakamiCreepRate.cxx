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

#include "neml2/models/solid_mechanics/LiuMurakamiCreepRate.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/cosh.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/linalg/eigh.h"

namespace neml2
{
register_NEML2_object(LiuMurakamiCreepRate);

OptionSet
LiuMurakamiCreepRate::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() =
      "Liu-Murakami creep rate model, for use with Liu-Murakami Damage Rate."
      "Creep rate given by \\f$ \\dot{\\varepsilon}_{c} = \\frac{3}{2} B \\sigma_{EQ}^{n-1} s_{ij}"
      "\\left[ \\frac{2(n+1)}{\\pi \\sqrt(1+3/n)} \\left) \\frac{\\sigma_{I}}{\\sigma_{EQ}} \\right)^{2}"
      "\\omega^{3/2} where: " ;  
      
  options.set_input("stress") = VariableName(STATE, "S");
  options.set("stress").doc() = "Stress";

  options.set_input("scalar_damage") = VariableName(STATE, "internal","w");
  options.set("scalar_damage").doc() = "Scalar Damage variable";

  options.set_output("creep_strain_rate") = VariableName(STATE, "internal", "Ec_rate");
  options.set("creep_strain_rate").doc() = "Rate of creep strain";

  options.set_parameter<TensorName<Scalar>>("scaling");
  options.set("scaling").doc() = "Scaling";

  options.set_parameter<TensorName<Scalar>>("exponent");
  options.set("exponent").doc() = "Exponent";


  return options;
}

LiuMurakamiCreepRate::LiuMurakamiCreepRate(const OptionSet & options)
  : Model(options),
    _w(declare_input_variable<Scalar>("scalar_damage")),
    _S(declare_input_variable<SR2>("stress")),
    _B(declare_parameter<Scalar>("B", "scaling", true)),
    _n(declare_parameter<Scalar>("n", "exponent", true)),
    _Ec_dot(declare_output_variable<SR2>("creep_strain_rate"))
{
}

void
LiuMurakamiCreepRate::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  // Get principal stresses
  const auto & [eigvals, eigvecs] = linalg::eigh(_S);
  const auto & s1 = Scalar(eigvals[-1]);

  // Get von mises equivalent stress
  const auto eps = machine_precision(_S.scalar_type());
  auto S = SR2(_S).dev();
  Scalar vm = std::sqrt(3.0 / 2.0) * S.norm(eps);
  double pi = 3.1415926535;

  if (out)
  {
    _Ec_dot =  SR2((3.0/2.0)*_B*pow(vm,_n-1)*S*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)));
  }

  if (dout_din)
  {
    if (_S.is_dependent()) //Right now not sure how to handle stress dependence
    //as damage stress derivative is a pain with s1 eigenvalue
      _Ec_dot.d(_S) = SR2::identity_map(S.options()) * Scalar((3.0/2.0)*_B*pow(vm,_n-1)*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)));

    if (_w.is_dependent())
      _Ec_dot.d(_w) = SR2((9*_B*(_n+1)*S*pow(s1,2.0)*sqrt(_w)*pow(vm,_n-3.0)*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)))/(2*pi*sqrt((_n+3/_n))));

    if (const auto * const B = nl_param("B"))
      _Ec_dot.d(*B) = SR2((3.0/2.0)*pow(vm,_n-1)*S*exp(((2*(_n+1))/(pi*sqrt(1+(3.0/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0)));

    if (const auto * const n = nl_param("n"))
      _Ec_dot.d(*n) = SR2(3.0*_B*S*pow(vm,_n-3)*exp(((2*(_n+1))/(pi*sqrt(1+(3/_n))))*pow(s1/vm,2.0)*pow(_w,3.0/2.0))*(pi*pow(_n,2.0)*pow((_n+3)/_n,3.0/2.0)*pow(vm,2.0)*log(vm) + (2*pow(_n,2.0)+9*_n +3)*pow(s1,2.0)*pow(_w,3.0/2.0) ) );

    // if (const auto * const q = nl_param("q"))
    //    _w_dot.d(*q) = (_A/pow(_q,2.0))*pow(sd,_p)*exp(_q*(_w-1))*(exp(_q)*(_q*_w -1)+_q + 1);
      
    // if (const auto * const alpha = nl_param("alpha"))
    //    _w_dot.d(*alpha) = (_A/_q)*_p*(exp(_q)-1)*exp(_q*(_w-1))*(vm-sd)*pow(sd,_p-1);

    //   }
}
} 
}// namespace neml2
