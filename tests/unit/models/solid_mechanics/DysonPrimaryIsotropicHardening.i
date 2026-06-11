[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/gamma_rate state/internal/k state/internal/s'
    input_Scalar_values = '0.1 50.0 10'
    output_Scalar_names = 'state/internal/k_rate'
    output_Scalar_values = '0.55'
  []
[]

[Models]
  [model]
    type = DysonPrimaryIsotropicHardening
    saturated_hardening = 100
    initial_hardening_rate = 110.0
  []
[]
