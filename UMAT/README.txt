This folder contains the material subroutines for ETFE nonlinear viscoelastic models.
- The file "template.cae" is an Abaqus cae file which already has the two nonlinear models included.
- The csv files contains the material parameters which are to be copied and pasted in the 
  Abaqus user material window in order to create such material. The number of State Variables is 163.
- The subroutines, to be linked to Abaqus when creating the job; in the template file, this was already performed.

The temperature must be input in Celsius, while all the units used elsewhere must be in SI. 
The analysis should be run with NONLINEAR GEOMETRY option active.
