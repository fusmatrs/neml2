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

#include "neml2/models/solid_mechanics/SinhPlasticFlowRate.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/heaviside.h"
#include "neml2/tensors/functions/abs.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/log.h"
#include "neml2/tensors/functions/sinh.h"
#include "neml2/tensors/functions/cosh.h"

namespace neml2
{
register_NEML2_object(SinhPlasticFlowRate);

OptionSet
SinhPlasticFlowRate::expected_options()
{
  OptionSet options = PlasticFlowRate::expected_options();
  options.doc() =
      "Hyperbolic sine viscous approximation of the consistent yield envelope "
      ", i.e. \\f$ \\dot{\\gamma} = \\alpha sinh\\left( \\frac{f}{\\beta} \\right) \\f$, "
      "where \\f$ f \\f$ is the yield function, \\f$ \\alpha \\f$ is the reference rate, and \\f$ "
      "\\beta \\f$ is the drag stress.";

  options.set_parameter<TensorName<Scalar>>("reference_rate");
  options.set("reference_rate").doc() = "Reference rate";

  options.set_parameter<TensorName<Scalar>>("drag_stress");
  options.set("drag_stress").doc() = "Drag stress";

  return options;
}

SinhPlasticFlowRate::SinhPlasticFlowRate(const OptionSet & options)
  : PlasticFlowRate(options),
    _alpha(declare_parameter<Scalar>("alpha", "reference_rate", /*allow_nonlinear=*/true)),
    _beta(declare_parameter<Scalar>("beta", "drag_stress", /*allow_nonlinear=*/true))
{
}

void
SinhPlasticFlowRate::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // Compute the Sinh approximation of the yield surface
  auto Hf = heaviside(Scalar(_f));
  auto f_abs = abs(Scalar(_f));
  auto gamma_dot_m = _alpha * sinh( f_abs / _beta );
  auto gamma_dot = gamma_dot_m * Hf;

  if (out)
    _gamma_dot = gamma_dot;

  if (dout_din)
  {
    auto dgamma_dot_df = (_alpha / _beta) * cosh( f_abs / _beta );

    if (_f.is_dependent())
      _gamma_dot.d(_f) = dgamma_dot_df;

    if (const auto * const alpha = nl_param("alpha"))
      _gamma_dot.d(*alpha) = sinh(f_abs / _beta);

    if (const auto * const beta = nl_param("beta"))
      _gamma_dot.d(*beta) = (-_alpha / (_beta * _beta)) * f_abs * cosh( f_abs / _beta );
  }
}
} // namespace neml2
