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

#include "neml2/models/solid_mechanics/damage/CavityEvolution.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/exp.h"

namespace neml2
{
register_NEML2_object(CavityEvolution);

OptionSet
CavityEvolution::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() = "Cavity evolution from Dyson";

  options.set_input("strain_rate") = VariableName(STATE, "internal","e_rate");
  options.set("strain_rate").doc() = "Strain Rate";

  options.set_parameter<TensorName<Scalar>>("growth_rate");
  options.set("growth_rate").doc() = "Growth Rate";
  options.set_parameter<TensorName<Scalar>>("failure_strain");
  options.set("failure_strain").doc() = "failure_strain";

  return options;
}

CavityEvolution::CavityEvolution(
    const OptionSet & options)
  : ScalarDamageRate(options),
    _e_rate(declare_input_variable<Scalar>("strain_rate")),
    _kn(declare_parameter<Scalar>("Kn", "growth_rate", true)),
    _efu(declare_parameter<Scalar>("efu", "failure_strain", true))

{
}

void
CavityEvolution::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  if (out)
  {
      _w_dot = (_kn/_efu)*_e_rate;
  }

  if (dout_din)
  {
    if (_e_rate.is_dependent())
      _w_dot.d(_e_rate) = (_kn/_efu);

    if (const auto * const Kn = nl_param("Kn"))
      _w_dot.d(*Kn) = (1/_efu)*_e_rate;

    if (const auto * const efu = nl_param("efu"))
      _w_dot.d(*efu) = -(_kn/pow(_efu,2))*_e_rate;
    
  }
}
} // namespace neml2
