[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/s state/w'
    input_Scalar_values = '100 0.5'
    output_Scalar_names = 'state/wdot'
    output_Scalar_values = '0.32000000000000006'
  []
[]

[Models]
  [model]
    type = LeckieHayhurstScalarDamageRate
    effective_stress = 'state/s'
    scalar_damage = 'state/w'
    scalar_damage_rate = 'state/wdot'
    stress_scaling = 500
    stress_exponent = 2
    damage_exponent = 3
  []
[]
