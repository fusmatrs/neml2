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

#include "neml2/models/solid_mechanics/LeckieHayhurstScalarDamageRate.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"

namespace neml2
{
register_NEML2_object(LeckieHayhurstScalarDamageRate);

OptionSet
LeckieHayhurstScalarDamageRate::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() = "Leckie Hayhurst Scalar Damage Rate, \\f$ \\dot{\\omega} = \\left( \\right)^\\zeta"
                  "\\left(1 - \\omega \\right)^-\\phi";

options.set_input("effective_stress") = VariableName(STATE,  "s");
  options.set("effective_stress").doc() = "Effective stress";

  options.set_parameter<TensorName<Scalar>>("stress_scaling");
  options.set("stress_scaling").doc() = "Stress Scaling";
  options.set_parameter<TensorName<Scalar>>("stress_exponent");
  options.set("stress_exponent").doc() = "Stress Exponent";
  options.set_parameter<TensorName<Scalar>>("damage_exponent");
  options.set("damage_exponent").doc() = "Damage Exponent";

  return options;
}

LeckieHayhurstScalarDamageRate::LeckieHayhurstScalarDamageRate(
    const OptionSet & options)
  : ScalarDamageRate(options),
    _s(declare_input_variable<Scalar>("effective_stress")),
    _A(declare_parameter<Scalar>("A", "stress_scaling", true)),
    _zeta(declare_parameter<Scalar>("zeta", "stress_exponent", true)),
    _phi(declare_parameter<Scalar>("phi", "damage_exponent", true))
{
}

void
LeckieHayhurstScalarDamageRate::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  const auto sp = pow(_s / _A,_zeta);
  const auto dp = pow(1-_w,-_phi);

  if (out)
    _w_dot = sp * dp;

  if (dout_din)
  {
    if (_s.is_dependent())
      _w_dot.d(_s) =  ( _zeta *  dp * sp ) / _s;

    if (_w.is_dependent())
      _w_dot.d(_w) = _phi * pow(1-_w,-_phi-1) * sp;

    if (const auto * const A = nl_param("A"))
      _w_dot.d(*A) = -( _zeta * dp * sp) / _A;

    if (const auto * const zeta = nl_param("zeta"))
      _w_dot.d(*zeta) = dp * sp * log(_s / _A);
    
    if (const auto * const phi = nl_param("phi"))
      _w_dot.d(*phi) = -dp * log(1-_w) * sp;
  }
}
} // namespace neml2
