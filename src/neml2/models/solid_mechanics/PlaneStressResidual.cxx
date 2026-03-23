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

#include "neml2/models/solid_mechanics/PlaneStressResidual.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/Scalar.h"

namespace neml2
{
register_NEML2_object(PlaneStressResidual);

OptionSet
PlaneStressResidual::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() =
      "Define the plane stress residual, i.e. s33=0";

  NonlinearSystem::enable_automatic_scaling(options);

  options.set_input("stress")= VariableName(STATE, "internal", "M");
  options.set("stress").doc() = "Current stress";

  options.set_input("e33_strain")= VariableName(STATE, "internal", "e33");
  options.set("e33_strain").doc() = "Out-of-plane strain";

  return options;
}

PlaneStressResidual::PlaneStressResidual(const OptionSet & options)
  : Model(options),
    _S(declare_input_variable<SR2>("stress")),
    _e33(declare_input_variable<Scalar>("e33_strain")),
    _r(declare_output_variable<Scalar>(_e33.name().remount(RESIDUAL)))
{
}


void
PlaneStressResidual::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto comps = at::split(_S, 1, -1);
  
  const auto & e33 = comps[2];

  if (out)
    _r = Scalar(e33.reshape(_e33.batch_sizes().concrete()),_e33.batch_sizes());

  if (dout_din)
  {
    _r.d(_S) = SR2::create({0, 0, 1, 0, 0, 0},_S.options());


    if (currently_solving_nonlinear_system())
      return;

    // _r.d(_sn) = -I;
    // _r.d(_t) = -_ds_dt;
    // _r.d(_tn) = _ds_dt;
  }
}

} // namespace neml2
