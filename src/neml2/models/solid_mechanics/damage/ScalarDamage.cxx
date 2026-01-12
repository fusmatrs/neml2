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

#include "neml2/models/solid_mechanics/damage/ScalarDamage.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/functions/pow.h"

namespace neml2
{
register_NEML2_object(ScalarDamage);

OptionSet
ScalarDamage::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Reduce effective stress by 1-damage ";

  options.set_input("effective_stress") = VariableName(STATE,  "s");
  options.set("effective_stress").doc() = "Effective stress";

  options.set_input("scalar_damage") = VariableName(STATE,  "w");
  options.set("scalar_damage").doc() = "Scalar Damage";

  options.set_output("damage_stress") = VariableName(STATE,  "d");
  options.set("damage_stress").doc() = "Damage Stress";

  return options;
}

ScalarDamage::ScalarDamage(const OptionSet & options)
  : Model(options),
    _s(declare_input_variable<Scalar>("effective_stress")),
    _w(declare_input_variable<Scalar>("scalar_damage")),
    _sd(declare_output_variable<Scalar>("damage_stress"))
{
}

void
ScalarDamage::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{

  auto wmod = Scalar::full(1.0) - _w; 

  if (out)
    _sd = _s / wmod;

  if (dout_din)
  {
    if (_s.is_dependent())
      _sd.d(_s)= 1.0 / wmod;
    
    if (_w.is_dependent())
      _sd.d(_w) = _s / pow(wmod,2.0);
  }
    
}
} // namespace neml2
