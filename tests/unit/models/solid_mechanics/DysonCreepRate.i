[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'forces/T state/internal/s state/internal/k state/internal/wp state/internal/wn state/internal/wd'
    input_Scalar_values = '900 10.0 0.5 0.0 0.0 0.0'
    input_SR2_names = 'state/S'
    input_SR2_values = 'S'
    output_SR2_names = 'state/internal/Ec_rate'
    output_SR2_values = 'er'
    check_derivatives = true
    #check_second_derivatives = true
  []
[]

[Tensors]
  [S]
    type = FillSR2
    values = '10.0 0.0 0.0'
  []
    [er]
    type = FillSR2
    values = '5210.256730803329 -2605.1283654016647 -2605.1283654016647'

  []
[]

[Models]
  [model]
    type = DysonCreepRate
    stress = 'state/S'
    temperature = 'forces/T'
    equivalent_stress='state/internal/s'
    isotropic_hardening = 'state/internal/k'
    dislocation_evolution = 'state/internal/wd'
    cavity_evolution = 'state/internal/wn'
    particle_evolution = 'state/internal/wp'
    activation_energy = 1
    rate_scaling = 1e4
    stress_scaling = 10.0
  []
[]
