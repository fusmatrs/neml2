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

#include "neml2/models/solid_mechanics/LinearScalarDamage.h"
#include "neml2/tensors/Scalar.h"

namespace neml2
{
register_NEML2_object(LinearScalarDamage);

OptionSet
LinearScalarDamage::expected_options()
{
  OptionSet options = ScalarDamage::expected_options();
  options.doc() += " following a linear relationship, i.e., \\f$ \\varomega = Kt \\f$ where "
                   "\\f$ K \\f$ is the damage modulus.";

  options.set<bool>("define_second_derivatives") = true;

  options.set_parameter<TensorName<Scalar>>("damage_modulus");
  options.set("damage_modulus").doc() = "Damage modulus";

  return options;
}

LinearScalarDamage::LinearScalarDamage(const OptionSet & options)
  : ScalarDamage(options),
    _K(declare_parameter<Scalar>("K", "damage_modulus"))
{
}

void
LinearScalarDamage::set_value(bool out, bool dout_din, bool d2out_din2)
{
  if (out)
    _w = _K + _w;

  if (dout_din)
    //if (_s.is_dependent())
    //  _w.d(_s) = 0;

  if (d2out_din2)
  {
    // zero
  }
}
} // namespace neml2
