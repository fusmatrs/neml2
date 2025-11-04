[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/fp'
    input_Scalar_values = '50'
    output_Scalar_names = 'state/internal/gamma_rate'
    output_Scalar_values = '0.0003395405572561501'
  []
[]

[Models]
  [model]
    type = SinhPlasticFlowRate
    reference_rate = 1e-3
    drag_stress = 150
  []
[]
