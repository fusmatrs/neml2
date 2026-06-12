[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/T state/internal/wn'
    input_Scalar_values = '100 0.1'
    output_Scalar_names = 'state/wndot'
    output_Scalar_values = '0.06553113671789303'
    check_derivatives = true
    #check_second_derivatives = true
  []
[]

[Models]
  [model]
    type = ParticleEvolution
    growth_rate=0.3
    temperature = 'state/T'
    activation_energy = 1
    scalar_damage = 'state/internal/wn'
    scalar_damage_rate = 'state/wndot'
    
  []
[]
