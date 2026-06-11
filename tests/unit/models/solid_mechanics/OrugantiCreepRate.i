[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/s state/T state/internal/ws state/internal/wp state/internal/k'
    input_Scalar_values = '100.0 900.0 0.0 0.0 4.0'
    input_SR2_names = 'state/S'
    input_SR2_values = 'S'
    output_SR2_names = 'state/internal/Ec_rate'
    output_SR2_values = 'er'
    check_derivatives = true
    derivative_abs_tol = 1#e0
    #check_second_derivatives = true
  []
[]

[Tensors]
  [S]
    type = FillSR2
    values = '100.0 0.0 0.0'
  []
    [er]
    type = FillSR2
    values = '4.682318932743124 -2.341159466371562 -2.341159466371562'
  []
[]

[Models]
  [model]
    type = OrugantiCreepRate
    stress = 'state/S'
    equivalent_stress='state/internal/s'
    temperature = 'state/T'
    mx_evolution = 'state/internal/wp'
    subgrain_evolution = 'state/internal/ws'
    isotropic_hardening = 'state/internal/k'
    rate_scaling = 1e-5
    stress_scaling = 0.45
    activation_energy = 1
  []
[]
