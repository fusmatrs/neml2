[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/s state/internal/w'
    input_Scalar_values = '100 0.5'
    output_Scalar_names = 'state/internal/fp'
    output_Scalar_values = '122.47448713915891'
  []
[]


[Models]

  [yield]
    type = YieldFunctionDamage
    yield_stress = 50
    
  []
  [model]
    type = ComposedModel
    models = 'yield'
  []
[]
