[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e_rate '
    input_Scalar_values = '1e-8'
    output_Scalar_names = 'state/internal/w'
    output_Scalar_values = '6e-9'
    check_derivatives = true
    #check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = CavityEvolution
    strain_rate = 'state/internal/e_rate'
    growth_rate=0.3
    failure_strain = 0.5
    
  []
[]
