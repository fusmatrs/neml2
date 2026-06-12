[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e_rate forces/T state/internal/wp'
    input_Scalar_values = '1e-8 900 0.1'
    output_Scalar_names = 'state/wpdot'
    output_Scalar_values = '0.020733229089222212'
    check_derivatives = true
    #check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = MXEvolution
    temperature = 'forces/T'
    scalar_damage = 'state/internal/wp'
    scalar_damage_rate = 'state/wpdot'
    growth_rate = 5
    initial_size =3.75
    activation_energy = 1
    
  []
[]
