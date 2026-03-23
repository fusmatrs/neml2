[Drivers]
  [unit]
    type = ModelUnitTest
    model = 'model'
    input_Scalar_names = 'state/internal/e33'
    input_Scalar_values = 'e33'
    #input_SR2_names = 'FORCES/E'
    #input_SR2_values = 'strains'#'E'
    output_SR2_names = 'state/internal/E_trial'
    output_SR2_values = 'out_strains'#'Et'
  []
[]

[Tensors]
  [E]
    type = FillSR2
    values = '10. 10. 10. 10. 10. 10.'
  []
  [Et]
    type = FillSR2
    values = '0. 0. 5. 0. 0. 0.'
  []
  [strains]
    type = LinspaceSR2
    start = 0
    end = E
    nstep = 5
  []
    [out_strains]
    type = LinspaceSR2
    start = 0
    end = Et
    nstep = 5
  []
  [e33]
  type = LinspaceScalar
  start=0
  end = 5
  nstep =5
  #batch_shape = (5,1)
  []

[]

[Models]
  [model]
    type = OutofPlaneStrain
    e33_strain = 'state/internal/e33'
  []
[]
