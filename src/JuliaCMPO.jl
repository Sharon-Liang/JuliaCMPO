module JuliaCMPO
__precompile__()

#=
### *Using Standard Libraries*
=#

using CUDA; CUDA.allowscalar(false)
using Parameters
using OMEinsum, LinearAlgebra
using LogExpFunctions
using Zygote, Optim, ChainRules
using Printf


#=
### *Includes And Exports* : *global.jl*
=#

#=
*Summary* :

Define some type aliases and string constants for the JuliaCMPO toolkit.

*Members* :

```text
Processor   -> Enumerated type for supported processors.
CPU, GPU    -> Values of Processor

#
__LIBNAME__ -> Name of this julia toolkit.
__VERSION__ -> Version of this julia toolkit.
__RELEASE__ -> Released date of this julia toolkit.
__AUTHORS__ -> Authors of this julia toolkit.
#
authors     -> Print the authors of JuliaCMPO to screen.
```
=#

#
include("global.jl")
#
export Processor, CPU, GPU 
#
export __LIBNAME__
export __VERSION__
export __RELEASE__
export __AUTHORS__
#
export authors



#=
### *Includes And Exports* : *types.jl*
=#

#=
*Summary* :

Define some dicts and structs, which are used to store the config
parameters or represent some essential data structures.

*Members* :

```text
OperatorType    -> Enumerated type for Fermionic/Bosonic correlators.
Fermi, Bose     -> Values of OperatorType
PauliMatrixName -> Enumerated type for Pauli matrices.
PX, PY, iPY, PZ, PPlus, PMinus  -> Values of PauliMatrixName.
#
CMPS          -> Data structure of cMPS local tensor.
CMPO          -> Data structure of cMPO local tensor.
```
=#

#
include("types.jl")
#
export OperatorType, Bose, Fermi
export PauliMatrixName, PX, PY, iPY, PZ, PPlus, PMinus
#
export CMPS
export CMPO



#=
### *Includes And Exports* : *utils.jl*
=#

#=
*Summary* :

Utilities

*Members* :

```text
symmetrize  -> Symmetrize a matrix.
diagm       -> Construct a square matrix of type `CuMatrix` form a `CuVector`. 
pauli       -> Generate Pauli matrices.
calc_ωn     -> Calculate Masubara frequency ωn.
#
eigensolver -> generate eigen values and vectors of a hermitian matrix.
#
logtrexp    -> logtrexp function.
#
optim_functions -> Generate gradient function for optimization.
#
solver_function -> Generate solver function 
cpu_solver      -> CPU solver function
gpu_solver      -> GPU solver function
```
=#


#
include("utils.jl")
#
export symmetrize
export pauli
export calc_ωn
#
export eigensolver
#
export logtrexp
#
export optim_functions
#
export solver_function
export cpu_solver
export gpu_solver



#=
### *Includes And Exports* : *math.jl*
=#

#=
*Summary* :

Mathematics for cMPO and cMPS local tensors.

*Members* :

```text
⊗           -> Multiplications between arrays in CMPS and CMPO data structures.
*           -> Multiplications between CMPS and CMPO structures.
#
log_overlap -> ln(⟨ψl|ψr⟩).
norm        -> √|⟨ψ|ψ⟩|.
normalize   -> normalize a cMPS.
logfidelity -> Calculate the logarithm of fidelity between two cMPS.
fidelity    -> Calculate fidelity between two cMPS.
project     -> Perform unitary transformation of a cMPS.
diagQ       -> Transform a cMPS to a gauge where Q is diagonalized.
#
transpose   -> Calculate the transpose of a cMPO.
adjoint     -> Calculate the adjoint of a cMPO.
#
==          -> Determine if two cMPS/cMPO are equal.
≈           -> Determine if two cMPS/cMPO are approximate.
ishermitian -> Determine if a cMPO is hermitian.
#
cat         -> cat two CMPO blocks


```
=#

#
include("math.jl")
#
export ⊗
#
export log_overlap
export logfidelity, fidelity
export project
export diagQ



#=
### *Includes And Exports* : *core.jl*
=#

#=
*Summary* :

Functions about initiating and compressing cMPS.

*Members* :

```text
MeraUpdateOptions    -> Keyword Arguments of mera_update function
_interpolate_isometry -> Interpolate between two unitrary matrices.
mera_update           -> Adaptive MERA update algorithm
#
CompressOptions -> Keyword Arguments of compress_cmps function
compress_cmps   -> To compress a cMPS.
#
init_cmps       -> initiate a cMPS
#


```
=#

#
include("core.jl")
#
export MeraUpdateOptions
export mera_update
#
export CompressOptions
export compress_cmps
#
export init_cmps



#=
### *Includes And Exports* : *rrule.jl*
=#

#=
*Summary* :

Rrules

*Members* :

```text
rrule(::typeof(logtrexp), t::Real, M::AbstractArray)
#


```
=#

#
include("rrule.jl")



#=
### *Includes And Exports* : *models.jl*
=#

#=
*Summary* :

CMPO Construction for various physical models.

*Members* :

```text
AbstractModel      -> Abstract type for physical model cMPO.
model              -> Construct model cMPO
#
TFIsingChain       -> Transcerse field Ising chain.
XXChain            -> 1D XX model.
XXZChain           -> 1D XXZ model.
J1J2Chain          -> 1D J1-J2 model.
#
TFIsingSquareHelical -> 2D transverse field ising model on a cylinder, square lattice, helical boundary condition.
XXSquareHelical      -> 2D XX model on a cylinder, square lattice, helical boundary condition.
XXZSquareHelical     -> 2D XXZ model on a cylinder, square lattice, helical boundary condition.
#


```
=#

#
include("models.jl")
#
export AbstractModel
export model
#
export TFIsingChain
export XXChain
export XXZChain
export J1J2Chain
#
export TFIsingSquareHelical
export XXSquareHelical
export XXZSquareHelical



#=
### *Includes And Exports* : *observables.jl*
=#

#=
*Summary* :

Physical observables.

*Members* :

```text
free_energy      -> Calculate free energy.
thermal_average  -> Calculate thermal average.
energy           -> Calculate energy density.
entropy          -> Calculate entropy.
specific_heat    -> Calculate specific heat.
#
correlation_2time -> Calculate C(τ) = -G(τ)
retarded_GF       -> Calculate G(z)
#


```
=#

#
include("observables.jl")
#
export free_energy
export thermal_average
export energy
export entropy
export specific_heat
#
export correlation_2time
export retarded_GF





end # module
