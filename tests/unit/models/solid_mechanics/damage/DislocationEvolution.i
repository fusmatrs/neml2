[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e_rate state/internal/wn'
    input_Scalar_values = '1 0.1'
    output_Scalar_names = 'state/wndot'
    output_Scalar_values = '0.243'
    check_derivatives = true
    #check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = DislocationEvolution
    strain_rate = 'state/internal/e_rate'
    constant_growth_rate=0.3
    #failure_strain = 0.5
    scalar_damage = 'state/internal/wn'
    scalar_damage_rate = 'state/wndot'
    
  []
[]
