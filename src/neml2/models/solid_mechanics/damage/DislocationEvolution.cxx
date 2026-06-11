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

#include "neml2/models/solid_mechanics/damage/DislocationEvolution.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/exp.h"

namespace neml2
{
register_NEML2_object(DislocationEvolution);

OptionSet
DislocationEvolution::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() =
      "Dislocation evolution from Dyson";

  options.set_input("strain_rate") = VariableName(STATE, "edot");
  options.set("strain_rate").doc() = "Equivalent Strain Rate";

  options.set_parameter<TensorName<Scalar>>("constant_growth_rate");
  options.set("constant_growth_rate").doc() = "Constant Growth Rate";

  return options;
}

DislocationEvolution::DislocationEvolution(const OptionSet & options)
  : ScalarDamageRate(options),
    _edot(declare_input_variable<Scalar>("strain_rate")),
    _C(declare_parameter<Scalar>("C", "constant_growth_rate", true))
{
}

void
DislocationEvolution::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  if (out)
  {
    _w_dot = (_edot*_C) * pow(1.0 - _w, 2.0);
  }

  if (dout_din)
  {

    if (_edot.is_dependent())
      _w_dot.d(_edot) =(1*_C) * pow(1.0 - _w, 2.0);

    if (_w.is_dependent())
      _w_dot.d(_w) = (-2.0*_edot * _C) * (1.0 - _w);

    if (const auto * const C = nl_param("C"))
      _w_dot.d(*C) = _edot * pow(1.0 - _w, 2.0);

  }
}
} // namespace neml2
