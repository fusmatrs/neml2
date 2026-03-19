[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e22'
    input_Scalar_values = '5.'
    input_SR2_names = 'FORCES/E'
    input_SR2_values = 'E'
    output_SR2_names = 'state/internal/E_trial'
    output_SR2_values = 'Et'
  []
[]

[Tensors]
  [E]
    type = FillSR2
    values = '10. 10. 10. 10. 10. 10.'
  []
  [Et]
    type = FillSR2
    values = '10. 10. 0. 0. 0. 10.'
  []
[]

[Models]
  [model]
    type = PlaneStress
    e33_strain = 'state/internal/e33'
    strain = 'FORCES/E'
  []
[]
