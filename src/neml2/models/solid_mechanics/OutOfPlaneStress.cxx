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

#include "neml2/models/solid_mechanics/OutOfPlaneStress.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/mandel_notation.h"
#include "neml2/tensors/functions/stack.h"

namespace neml2
{
register_NEML2_object(OutOfPlaneStress);

OptionSet
OutOfPlaneStress::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Spit out the out of plane stress s33";

  options.set_input("stress") = VariableName(STATE, "S");
  options.set("stress").doc() = "Stress";

  options.set_output("stress_33") = VariableName(STATE, "internal", "E_trial");
  options.set("stress_33").doc() = "Out of plane stress";

  return options;
}

OutOfPlaneStress::OutOfPlaneStress(const OptionSet & options)
  : Model(options),
    _S(declare_input_variable<SR2>("stress")),
    _s33(declare_output_variable<Scalar>("stress_33"))
{
}

void
OutOfPlaneStress::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  if (out)
  {
      // Copy strain
      const auto comps = at::split(_S, 1, -1);
      _s33 = Scalar(comps[2].reshape(_s33.batch_sizes().concrete()),_s33.batch_sizes());
 
}

  if (dout_din)
  {
    _s33.d(_S) = SR2::create({0, 0, 1, 0, 0, 0},_S.options());
        
  }
}
} // namespace neml2
