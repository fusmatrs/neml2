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

#include "neml2/models/solid_mechanics/OutofPlaneStrainResidual.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/mandel_notation.h"
#include "neml2/tensors/functions/stack.h"

namespace neml2
{
register_NEML2_object(OutofPlaneStrainResidual);

OptionSet
OutofPlaneStrainResidual::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Assuming ischoric flow, generate out of plane strain";

  NonlinearSystem::enable_automatic_scaling(options);

  options.set_input("e33_strain") = VariableName(STATE, "internal" "e33");
  options.set("e33_strain").doc() = "Out-of-plane strain";

  options.set_input("elastic_strain") = VariableName(STATE, "internal", "Ee");
  options.set("elastic_strain").doc() = "Elastic Strain";

  options.set_input("plastic_strain") = VariableName(STATE, "internal", "Ep");
  options.set("plastic_strain").doc() = "Plastic Strain";

  options.set_parameter<TensorName<Scalar>>("poisson_ratio");
  options.set("poisson_ratio").doc() = "Poissons Ratio";

  return options;
}

OutofPlaneStrainResidual::OutofPlaneStrainResidual(const OptionSet & options)
  : Model(options),
    _e33(declare_input_variable<Scalar>("e33_strain")),
    _Ee(declare_input_variable<SR2>("elastic_strain")),
    _Ep(declare_input_variable<SR2>("plastic_strain")),
    _nu(declare_parameter<Scalar>("nu", "poisson_ratio", true)),
    _r(declare_output_variable<Scalar>(_e33.name().remount(RESIDUAL)))
{
}

void
OutofPlaneStrainResidual::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  const auto zero = Scalar::zeros_like(_e33);
  const auto one = Scalar::ones_like(_e33);
  const auto msk = SR2::fill(one,one,zero,zero,zero,zero);
  const auto el = SR2(_Ee).inner(msk);
  const auto pl = SR2(_Ep).inner(msk);
  if (out)
  {
    _r = -(_nu/(1-_nu))*el - pl - _e33;

}

  if (dout_din)
  {
    if (_e33.is_dependent())      
      _r.d(_e33) = -one;

    if (_Ee.is_dependent())
      _r.d(_Ee) = -(_nu/(1-_nu))*msk;

    if (_Ep.is_dependent())
      _r.d(_Ep) = -msk;

    if (const auto * const nu = nl_param("nu"))
      _r.d(*nu) = -1/(1-_nu)*el;

        
  }

}
} // namespace neml2
