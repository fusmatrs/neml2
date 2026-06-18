[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e33 '
    input_Scalar_values = 'e33'
    input_SR2_names = 'state/internal/Ee state/internal/Ep'
    input_SR2_values = 'Ee Ep'
    output_Scalar_names = 'residual/internal/e33'
    output_Scalar_values = 'res'#'Et'
  []
[]

[Tensors]
  [E]
    type = FillSR2
    values = '2. 2. 2. 0. 0. 0.'
  []
  [Ee]
    type = FillSR2
    values = '2. 2. 2. 0. 0. 0.'
  []

  [Ep]
  type = FillSR2
  values = '1. 1. 1. 0. 0. 0.'
  []
  
  [e33]
    type = Scalar
    values=1
  []

  [res]
    type = Scalar
    values=-4.3333333333
  []

[]

[Models]
  [model]
    type = OutofPlaneStrain
    e33_strain = 'state/internal/e33'
    elastic_strain = 'state/internal/Ee'
     plastic_strain = 'state/internal/Ep'
     poisson_ratio = 0.25
  []
[]
