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

#include "neml2/models/solid_mechanics/damage/SubgrainEvolution.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/exp.h"

namespace neml2
{
register_NEML2_object(SubgrainEvolution);

OptionSet
SubgrainEvolution::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() =
      "Subgrain evolution from Oruganti et al, \\f$ \\dot{\\omega} = \\left( \\right)^\\zeta"
      "\\left(1 - \\omega \\right)^-\\phi";

  options.set_input("temperature") = VariableName(STATE, "T");
  options.set("temperature").doc() = "Temperature in K";

  options.set_input("strain_rate") = VariableName(STATE, "edot");
  options.set("strain_rate").doc() = "Equivalent Strain Rate";

  options.set_parameter<TensorName<Scalar>>("constant_growth_rate");
  options.set("constant_growth_rate").doc() = "Constant Growth Rate";
  options.set_parameter<TensorName<Scalar>>("exponential_growth_rate");
  options.set("exponential_growth_rate").doc() = "Exponential Growth Rate";
  options.set_parameter<TensorName<Scalar>>("initial_size");
  options.set("initial_size").doc() = "Initial Size";
  options.set_parameter<TensorName<Scalar>>("activation_energy");
  options.set("activation_energy").doc() = "Activation Energy";

  return options;
}

SubgrainEvolution::SubgrainEvolution(const OptionSet & options)
  : ScalarDamageRate(options),
    _T(declare_input_variable<Scalar>("temperature")),
    _edot(declare_input_variable<Scalar>("strain_rate")),
    _ks1(declare_parameter<Scalar>("Ks1", "constant_growth_rate", true)),
    _ks2(declare_parameter<Scalar>("Ks2", "exponential_growth_rate", true)),
    _si(declare_parameter<Scalar>("Si", "initial_size", true)),
    _qs(declare_parameter<Scalar>("Qp", "activation_energy", true))
{
}

void
SubgrainEvolution::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // Gas Constant
  const auto R = 8.31446261815324; // JK-1mol-1
  auto et = exp(-_qs / (R * _T()));

  if (out)
  {
    _w_dot = (_edot / _si) * (_ks1 + _ks2 * et) * pow(1.0 - _w, 2.0);
  }

  if (dout_din)
  {
    if (_T.is_dependent())
      _w_dot.d(_T) = (_edot / (_si*R*pow(_T(),2.0))) * _ks2 * et * _qs * pow(1.0 - _w, 2.0);

    if (_edot.is_dependent())
      _w_dot.d(_edot) =(1.0/ _si) * (_ks1 + _ks2 * et) * pow(1.0 - _w, 2.0);

    if (_w.is_dependent())
      _w_dot.d(_w) = (-2.0*_edot / _si) * (_ks1 + _ks2 * et) * (1.0 - _w);

    if (const auto * const Ks1 = nl_param("Ks1"))
      _w_dot.d(*Ks1) = (_edot / _si) * pow(1.0 - _w, 2.0);

    if (const auto * const Ks2 = nl_param("Ks2"))
      _w_dot.d(*Ks2) = (-_edot / _si) * et * pow(1.0 - _w, 2.0);

    if (const auto * const Si = nl_param("Si"))
      _w_dot.d(*Si) = (-_edot / pow(_si,2.0)) * (_ks1 + _ks2 * et) * pow(1.0 - _w, 2.0);

    if (const auto * const Qs = nl_param("Qs"))
      _w_dot.d(*Qs) = (_edot / (_si*R*_T())) * _ks2 * et * pow(1.0 - _w, 2.0);
  }
}
} // namespace neml2
