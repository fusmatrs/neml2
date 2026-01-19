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

#include "neml2/models/solid_mechanics/damage/MXEvolution.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/exp.h"

namespace neml2
{
register_NEML2_object(MXEvolution);

OptionSet
MXEvolution::expected_options()
{
  OptionSet options = ScalarDamageRate::expected_options();
  options.doc() = "MX evolution from Oruganti et al (2011), \\f$ \\dot{D}_{P} = \\frac{K_{P}}{P_{i}^{3}}"
      "e^{-Q_{P}/RT}\\left(1- D_{P} \\right)^{4} \\f$";


options.set_input("temperature") = VariableName(STATE,  "T");
  options.set("temperature").doc() = "Temperature in K";

  options.set_parameter<TensorName<Scalar>>("growth_rate");
  options.set("growth_rate").doc() = "Growth Rate";
  options.set_parameter<TensorName<Scalar>>("initial_size");
  options.set("initial_size").doc() = "Initial Size";
  options.set_parameter<TensorName<Scalar>>("activation_energy");
  options.set("activation_energy").doc() = "Activation Energy";

  return options;
}

MXEvolution::MXEvolution(
    const OptionSet & options)
  : ScalarDamageRate(options),
    _T(declare_input_variable<Scalar>("temperature")),
    _kp(declare_parameter<Scalar>("Kp", "growth_rate", true)),
    _pi(declare_parameter<Scalar>("Pi", "initial_size", true)),
    _qp(declare_parameter<Scalar>("Qp", "activation_energy", true))
{
}

void
MXEvolution::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // Gas Constant
  const auto R = 8.31446261815324; //JK-1mol-1
  auto et = exp(-_qp/(R*_T));

  if (out)
  {
      _w_dot = (_kp/(3.0*pow(_pi,3.0)))*et*pow(1.0-_w,4.0);
  }

  if (dout_din)
  {
    if (_T.is_dependent())
      _w_dot.d(_T) =  ((_kp*_qp)/(3.0*pow(_pi,3.0)*R*pow(_T,2.0)))*et*pow(1.0-_w,4.0);

    if (_w.is_dependent())
      _w_dot.d(_w) = ((-4.0*_kp)/(3.0*pow(_pi,3.0)))*et*pow(1.0-_w,3.0);

    if (const auto * const Kp = nl_param("Kp"))
      _w_dot.d(*Kp) = (1.0/(3.0*pow(_pi,3.0)))*et*pow(1.0-_w,4.0);

    if (const auto * const Pi = nl_param("Pi"))
      _w_dot.d(*Pi) = -3.0*(_kp/(3.0*pow(_pi,4.0)))*et*pow(1.0-_w,4.0);
    
    if (const auto * const Qp = nl_param("Qp"))
      _w_dot.d(*Qp) = (-_kp/(3.0*pow(_pi,3.0)*R*_T))*et*pow(1.0-_w,4.0);
  }
}
} // namespace neml2
