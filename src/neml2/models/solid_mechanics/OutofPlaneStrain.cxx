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

#include "neml2/models/solid_mechanics/OutofPlaneStrain.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/mandel_notation.h"
#include "neml2/tensors/functions/stack.h"

namespace neml2
{
register_NEML2_object(OutofPlaneStrain);

OptionSet
OutofPlaneStrain::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Produce an SR2 correction of 0,0,e33, 0, 0 ,0 .";

  options.set_input("e33_strain") = VariableName(STATE, "internal" "e33");
  options.set("e33_strain").doc() = "Out-of-plane strain";

  options.set_output("trial_strain") = VariableName(STATE, "internal", "E_trial");
  options.set("trial_strain").doc() = "Trial Strain";

  return options;
}

OutofPlaneStrain::OutofPlaneStrain(const OptionSet & options)
  : Model(options),
    _e33(declare_input_variable<Scalar>("e33_strain")),
    _Et(declare_output_variable<SR2>("trial_strain"))
{
}

void
OutofPlaneStrain::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto zero = Scalar::zeros_like(_e33);
  const auto one = Scalar::ones_like(_e33);
  if (out)
  {
    
   _Et = SR2::fill(zero,zero,_e33,zero,zero,zero);
  //_Et = SR2(base_stack({e11, e22, e33, zero, zero, e12}, -1));


}

  if (dout_din)
  {
    if (_e33.is_dependent())      
      _Et.d(_e33) = SR2::fill(zero,zero,one,zero,zero,zero);
        
  }

}
} // namespace neml2
