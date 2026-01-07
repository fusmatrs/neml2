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

#include "neml2/models/solid_mechanics/LiuMurakamiScalarDamageRate.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/heaviside.h"
#include "neml2/tensors/functions/linalg/eigh.h"

namespace neml2
{
register_NEML2_object(LiuMurakamiScalarDamageRate);

OptionSet
LiuMurakamiScalarDamageRate::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() = "Liu Murakami Scalar Damage Rate, \\f$ \\dot{\\omega} = \\frac{A(1-e^{-q})}{q} "
                  "\\sigma_{D}^{p} e^{q} \\omega} \\f$";

  options.set_input("stress") = VariableName(STATE, "S");
  options.set("stress").doc() = "Stress";

  options.set_parameter<TensorName<Scalar>>("scaling");
  options.set("scaling").doc() = "Scaling";
  options.set_parameter<TensorName<Scalar>>("exponent");
  options.set("exponent").doc() = "Exponent";
  options.set_parameter<TensorName<Scalar>>("stress_exponent");
  options.set("stress_exponent").doc() = "Stress Exponent";
  options.set_parameter<TensorName<Scalar>>("stress_weight");
  options.set("stress_weight").doc() = "Stress Weight";

  return options;
}

LiuMurakamiScalarDamageRate::LiuMurakamiScalarDamageRate(const OptionSet & options)
  : ScalarDamageRate(options),
    _S(declare_input_variable<SR2>("stress")),
    _A(declare_parameter<Scalar>("A", "scaling", true)),
    _q(declare_parameter<Scalar>("q", "exponent", true)),
    _p(declare_parameter<Scalar>("p", "stress_exponent", true)),
    _alpha(declare_parameter<Scalar>("alpha", "stress_weight", true))
{
}

void
LiuMurakamiScalarDamageRate::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  // Get principal stresses
  const auto & [eigvals, eigvecs] = linalg::eigh(_S);
  const auto & s1 = Scalar(eigvals[-1]);

  // Get von mises equivalent stress
  const auto eps = machine_precision(_S.scalar_type());
  auto S = SR2(_S).dev();
  Scalar vm = std::sqrt(3.0 / 2.0) * S.norm(eps);

  // Get the 'Damage stress'
  Scalar sd = _alpha*vm + (1-_alpha)*s1;

  if (out)
  {
    _w_dot =  (_A/_q)*(1-exp(-_q))*pow(sd,_p)*exp(_q*_w);
  }

  if (dout_din)
  {
    // if (_s.is_dependent()) Right now not sure how to handle stress dependence
    // as damage stress derivative is a pain with s1 eigenvalue
    //   _w_dot.d(_s) = (_zeta * dp * sp) / _s;

    if (_w.is_dependent())
       _w_dot.d(_w) = _A*(exp(_q)-1)*exp(_q*(_w-1))*pow(sd,_p);

    if (const auto * const A = nl_param("A"))
       _w_dot.d(*A) = (1/_q)*(1-exp(-_q))*pow(sd,_p)*exp(_q*_w);

    if (const auto * const p = nl_param("p"))
       _w_dot.d(*p) = (_A/_q)*(exp(_q)-1)*pow(sd,_p)*exp(_q*(_w-1))*log(sd);

    if (const auto * const q = nl_param("q"))
       _w_dot.d(*q) = (_A/pow(_q,2.0))*pow(sd,_p)*exp(_q*(_w-1))*(exp(_q)*(_q*_w -1)+_q + 1);
      
    if (const auto * const alpha = nl_param("alpha"))
       _w_dot.d(*alpha) = (_A/_q)*_p*(exp(_q)-1)*exp(_q*(_w-1))*(vm-sd)*pow(sd,_p-1);

      }
}
} // namespace neml2
