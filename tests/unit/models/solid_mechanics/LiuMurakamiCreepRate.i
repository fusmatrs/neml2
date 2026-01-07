[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/w'
    input_Scalar_values = '0.1'
    input_SR2_names = 'state/S'
    input_SR2_values = 'S'
    output_SR2_names = 'state/internal/Ec_rate'
    output_SR2_values = 'Ecr'
  []
[]

[Tensors]
  [S]
    type=FillSR2
    values='100.0 50.0 50.0 0 0 0'
  []
  [Ecr]
    type=FillSR2
    values='54.192885977071334 -27.09633731239744 -27.09633731239744 0 0 0'
  []
[]

[Models]
  [model]
    type = LiuMurakamiCreepRate
    stress = 'state/S'
    scalar_damage = 'state/internal/w'
    scaling = 1
    exponent = 1
  []

[]
