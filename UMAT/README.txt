This folder contains the material subroutines for ETFE nonlinear viscoelastic models.
- Three material models are provided implemented in UMATs: Nonlinear Viscoelastic Isotropic model (NLVEEYISO), Nonlinear Viscoelastic Orthotropic 
  model (NLVEEYORTHO) and Nonlinear Viscoelastic Isotropic with Strain Rate Temperature Dependent Yield Criterion model (NLVEEYISO_YIELD).
- The file "template.cae" is an Abaqus cae file which already has the two nonlinear models included (isotropic and orthotropic without yield criterion).
- The csv files contains the material parameters which are to be copied and pasted in the 
  Abaqus user material window in order to create such material. The number of State Variables is 163 without yield criterion and 167 with yield criterion.
- The subroutines, to be linked to Abaqus when creating the job; in the template file, this was already performed.
- The folder "Example_Isotropic_Yield" has an example in Abaqus of a rectangular geometry uniaxially loaded modelled with NLVEEYISO_YIELD. 
  The files provided inside the folder are: input file and odb file with results

The temperature must be input in Celsius, while all the units used elsewhere must be in SI. 
The analysis should be run with NONLINEAR GEOMETRY option active.