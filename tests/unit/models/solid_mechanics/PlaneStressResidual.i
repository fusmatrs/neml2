[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_SR2_names = 'state/internal/M'
    input_SR2_values = 'stress'
    input_Scalar_names = 'state/internal/e33'
    input_Scalar_values ='e33'
    output_Scalar_names = 'residual/internal/e33'
    output_Scalar_values = 'out_stress'

  []
[]

[Tensors]
  [E]
    type = FillSR2
    values = '10. 10. 10. 10. 10. 10.'
  []
  [Et]
    type = FillSR2
    values = '10. 10. 0. 0. 0. 10.'
  []
  [stress]
    type = LinspaceSR2
    start = 0
    end = E
    nstep = 5
  []
    [out_stress]
    type = LinspaceScalar
    start = 0
    end = 10
    nstep = 5
  []
[e33]
  type = Scalar
  values=1
  batch_shape = (5)
[]
[]

[Models]
  [model]
    type = PlaneStressResidual
    stress = 'state/internal/M'
    
  []
[]
