############################################
# QuantumOutputDefinitions

"""
	abstract type State

Abstract base type for all quantum state representations.
"""
abstract type State end

""" 
	struct QuantumOutput{Z, I, S <: State}

Represents the complete results of a quantum program simulation, capturing the evolution of quantum states across multiple operations and shots.

This struct is **compile-time specialized** on the type parameters `Z` and `I` (that encode the qubit indexing settings) and on type parameter `S <: State` (that forces the struct to hold type-stable quantum state objects), ensuring optimized and branch-free code generation.

# Type parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.
- `S <: State`: The concrete type of the quantum states being stored, e.g., `VectorState` or `DensityState`.

# Fields
- `settings::Settings{Z, I}` — Qubit indexing settings.
- `numberOfQubits::Int` — Total number of qubits.
- `numberOfStates::Int` — Number of independent quantum states (typically 1 initial state + 1 for each operation).
- `numberOfShots::Int` — Number of repeated measurements (shots) per state.
- `output::Matrix{S}` — A matrix of size `(numberOfStates, numberOfShots)` containing the quantum state objects.

# Notes
Each element `output[i, j]` is a subtype of `State`, such as `VectorState` or `DensityState`.

# See also
- [`Settings`](@ref)
- [`State`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)

"""
struct QuantumOutput{Z, I, S <: State}
	settings::Settings{Z, I}
	numberOfQubits::Int
	numberOfStates::Int
	numberOfShots::Int
	output::Matrix{S}
end

"""
    struct Measured

Represents the result of a quantum measurement operation.

Let ``M`` be the number of measured qubits.  
The measurement result includes probabilities for all ``2^M`` possible outcomes, the selected outcome, and the corresponding binary labels.

# Fields
- `probability::Vector{Float64}` — A vector of length ``2^M`` containing the probabilities of each possible outcome.
- `outcome::Union{Int64, Nothing}` — The observed measurement result as an integer in the range ``0:(2^M - 1)``, or `nothing` if the result is not revealed.
- `labels::Vector{String}` — A vector of binary string labels corresponding to each outcome (e.g., `"00"`, `"01"`, etc.).

# Example
Extract the measured result from a quantum state of a single qubit in a perfect superposition:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> measureGate!(qc, [1])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> getStateShot(qo, 3, 1).measured
Measured([0.5, 0.5], 1, ["0", "1"])
```
"""
struct Measured
	probability::Vector{Float64}
	outcome::Union{Int64, Nothing}
	labels::Vector{String}
end

""" 
    struct VectorState <: State

Represents a pure quantum vector state in vector (ket) form: ``|q⟩``.

# Fields
- `q::Matrix{ComplexF64}` — Column vector representing the quantum vector state,  
  where ``|q⟩ = q₀|00...0⟩ + q₁|00...1⟩ + ... + q_{N-1}|11...1⟩`` with length ``N = 2^M`` with ``M`` equal to the number of qubits.
- `measured::Union{Measured, Nothing}` — The measurement result, if a measurement was performed.

# Notes
`VectorState` is broadcastable.

# See also
- [`Measured`](@ref)
- [`State`](@ref)

# Example
Extract the quantum vector state after applying a Hadamard gate to a single-qubit zero quantum vector state:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> getStateShot(qo, 2, 1).q
2×1 Matrix{ComplexF64}:
 0.7071067811865476 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
mutable struct VectorState <: State
	q::ComplexMatrix
	measured::Union{Measured, Nothing}
	VectorState(q::ComplexMatrix; measured = nothing) = new(q, measured)
end

VectorState(q::AbstractVector{<:Number}; measured = nothing) = VectorState(reshape(ComplexF64.(q), :, 1), measured=measured)

function VectorState(q::AbstractMatrix{<:Number}; measured = nothing)
	if size(q)[2] == 1
		return VectorState(ComplexF64.(q), measured = measured)
	else
		throw("The number of columns of q is not equal to one.")
	end
end

Base.broadcastable(s::VectorState) = Ref(s)

"""
    struct DensityState <: State

Represents a mixed quantum density state in matrix form: ``ρ``.

# Fields
- `rho::Matrix{ComplexF64}` — Quantum density state ``ρ = ∑ ρ_{ij} |i⟩⟨j|`` of size ``N × N``, where ``N = 2^M`` with ``M`` equal to the number of qubits.
- `measured::Union{Measured, Nothing}` — The measurement result, if a measurement was performed.

# Notes
- `DensityState` is broadcastable.

# See also
- [`Measured`](@ref)
- [`State`](@ref)

# Example
Extract the quantum density state after applying a Hadamard gate to a single-qubit zero pure quantum density state:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(densityFast, [[1.00, [1. 0.;]]])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> getStateShot(qo, 2, 1).rho
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.5+0.0im
 0.5+0.0im  0.5+0.0im
```
"""
mutable struct DensityState <: State
	rho::ComplexMatrix
	measured::Union{Measured, Nothing}
	DensityState(rho::ComplexMatrix; measured = nothing) = new(rho, measured)
end

DensityState(rho::AbstractMatrix{<:Number}; measured = nothing) = DensityState(ComplexF64.(rho), measured = measured)

Base.broadcastable(s::DensityState) = Ref(s)

""" 
    getState(
		quantumOutput::QuantumOutput, 
		stateId::Int
	) -> Vector{State}

Extracts all measurement shots for a given quantum state from the `quantumOutput`.

# Arguments
- `quantumOutput::QuantumOutput` — The result of running a quantum program.
- `stateId::Int` — Index of the quantum state to extract (1-based or 0-based depending on `settings`).

# Returns
- `Vector{State}` — A vector of length `numberOfShots` containing the state for each shot.

# See also
- [`State`](@ref)
- [`QuantumOutput`](@ref)
- [`Settings`](@ref)

# Example
Extract state 3 from a quantum output with 5 shots:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> measureGate!(qc, [1])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 5)
julia> getState(qo, 3)
5-element Vector{VectorState}:
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
```
"""
function getState(quantumOutput::QuantumOutput, stateId::Int)
	return quantumOutput.output[convertToOneBasedNumbering(quantumOutput.settings, stateId), :] # : for the numberOfShots
end

""" 
    getShot(
		quantumOutput::QuantumOutput, 
		shotId::Int
	) -> Vector{State}

Extracts all quantum states for a specific measurement shot from the `quantumOutput`.

# Arguments
- `quantumOutput::QuantumOutput` — The result of executing a quantum program.
- `shotId::Int` — Index of the shot to extract (1-based or 0-based depending on `settings`).

# Returns
- `Vector{State}` — A vector of length `numberOfStates` representing the system's states for the selected shot.

# See also
- [`State`](@ref)
- [`QuantumOutput`](@ref)
- [`Settings`](@ref)

# Example
Extract shot 4 from a quantum output with 3 states:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> measureGate!(qc, [1])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 5)
julia> getShot(qo, 4)
3-element Vector{VectorState}:
 VectorState(ComplexF64[1.0 + 0.0im; 0.0 + 0.0im;;], nothing)
 VectorState(ComplexF64[0.7071067811865476 + 0.0im; 0.7071067811865475 + 0.0im;;], nothing)
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
```
"""
function getShot(quantumOutput::QuantumOutput, shotId::Int)
	return quantumOutput.output[:,convertToOneBasedNumbering(quantumOutput.settings, shotId)] # : for the numberOfStates
end

""" 
    getStateShot(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		shotId::Int
	) -> State

Extracts the quantum state corresponding to a specific state index and shot index from the `quantumOutput`.

# Arguments
- `quantumOutput::QuantumOutput` — The result of executing a quantum program.
- `stateId::Int` — Index of the desired state (1-based or 0-based depending on `settings`).
- `shotId::Int` — Index of the desired shot (1-based or 0-based depending on `settings`).

# Returns
- `State` — The quantum state at the specified state and shot combination.

# See also
- [`State`](@ref)
- [`QuantumOutput`](@ref)
- [`Settings`](@ref)

# Example
Extract state 3 and shot 4 from a quantum output:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> measureGate!(qc, [1])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 5)
julia> getStateShot(qo, 3, 4)
VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
```
"""
function getStateShot(quantumOutput::QuantumOutput, stateId::Int, shotId::Int)
	return quantumOutput.output[convertToOneBasedNumbering(quantumOutput.settings, stateId), convertToOneBasedNumbering(quantumOutput.settings, shotId)]
end

""" 
    convertVectorStateToDensityState(
		vectorState::VectorState
	) -> DensityState

Converts a pure quantum vector state ``|q⟩`` into a quantum density state ``ρ = |q⟩⟨q|``.

# Arguments
- `vectorState::VectorState` — A quantum vector state ``|q⟩`` represented by a normalized complex vector.

# Returns
- `DensityState` — The corresponding quantum density state representation ``ρ`` of the pure quantum vector state.

# See also
- [`VectorState`](@ref)
- [`DensityState`](@ref)

# Example
Convert a vector state, holding a single qubit in superposition, into a density state:
```julia-repl
julia> vectorState = VectorState(reshape([1.; 1.]/sqrt(2).+0im, 2, 1))
julia> vectorState.q
2×1 Matrix{ComplexF64}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
julia> densityState = convertVectorStateToDensityState(vectorState)
DensityState(ComplexF64[0.4999999999999999 + 0.0im 0.4999999999999999 + 0.0im; 0.4999999999999999 - 0.0im 0.4999999999999999 + 0.0im], nothing)
julia> densityState.rho
2×2 Matrix{ComplexF64}:
 0.5+0.0im  0.5+0.0im
 0.5-0.0im  0.5+0.0im
```
"""
function convertVectorStateToDensityState(vectorState::VectorState)
	return DensityState(vectorState.q * vectorState.q')
end
