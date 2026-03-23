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

#include "neml2/models/solid_mechanics/PlaneStress.h"
#include "neml2/tensors/SR2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/mandel_notation.h"
#include "neml2/tensors/functions/stack.h"

namespace neml2
{
register_NEML2_object(PlaneStress);

OptionSet
PlaneStress::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Modify strains to enforce plane stress, requires PlaneStressResidual too.";

  options.set_input("strain") = VariableName(FORCES, "E");
  options.set("strain").doc() = "Strain";

  options.set_input("e33_strain") = VariableName(STATE, "internal" "e33");
  options.set("e33_strain").doc() = "Out-of-plane strain";

  options.set_output("trial_strain") = VariableName(STATE, "internal", "E_trial");
  options.set("trial_strain").doc() = "Trial Strain";

  return options;
}

PlaneStress::PlaneStress(const OptionSet & options)
  : Model(options),
    _E(declare_input_variable<SR2>("strain")),
    _e33(declare_input_variable<Scalar>("e33_strain")),
    _Et(declare_output_variable<SR2>("trial_strain"))
{
}

void
PlaneStress::set_value(bool out, bool dout_din, bool d2out_din2)
{

  if (out)
  {
      // Copy strain
      const auto comps = at::split(_E, 1, -1);
  //     const auto & e11 = comps[0];
  //     const auto & e22 = comps[1];
  //     const auto & e12 = comps[5];
  //     const auto zero = at::zeros_like(e11);
  //     const auto e33 = at::ones_like(e11) * _e33 ;

  // _Et.base_index_put_()

//   // set e22 as the trial and set _e13 and _e23 to zero. Assumes isotropy
  const auto & e11 = Scalar(comps[0].reshape(_e33.batch_sizes().concrete()),_e33.batch_sizes());
  const auto & e22 = Scalar(comps[1].reshape(_e33.batch_sizes().concrete()),_e33.batch_sizes());
  const auto & e12 = Scalar(comps[5].reshape(_e33.batch_sizes().concrete()),_e33.batch_sizes());
///mandel_factor(5);

  
  const auto zero = Scalar::zeros_like(e11);
  const auto e33 = Scalar::ones_like(e11) * _e33 ;
  //_Et = SR2::fill(Scalar(e00),Scalar(e11),_e22,zero,zero,Scalar(e01));
  _Et = SR2(base_stack({e11, e22, e33, zero, zero, e12}, -1));
  // Try assembling as torch objects. A la SR2 currently gives incorrect batch dims
  // const auto E = at::cat({e11,
  //                         e22,
  //                         e33,
  //                         zero,
  //                         zero,
  //                         e12},
  //                         -1);
  // _Et = SR2(E,_E.batch_sizes());
  // Alternative approach to get ezz as a variable too.

}

  if (dout_din)
  {
    if (_E.is_dependent())   
      _Et.d(_E) =  SSR4::create({{1, 0, 0, 0, 0, 0},
                      {0, 1, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0, 0},
                      {0, 0, 0, 0, 0, 1}},
                      _E.options());
    
    if (_e33.is_dependent())
        _Et.d(_e33) = SR2::create({0, 0, 1, 0, 0, 0},_E.options());
        
  }
    if (d2out_din2)
  {
    // zero
  }
}
} // namespace neml2
