############################################
# QuantumCircuitDefinitions

"""
	Qubits = Vector{Int}

Alias for a vector of integers representing qubit indices in a quantum circuit.
"""
Qubits = Vector{Int}

"""
	abstract type Gate

Abstract type representing any quantum gate. Concrete subtypes include `UnitaryGate`, `MeasureGate`, and `QuantumChannelGate`.

# See also
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)
"""
abstract type Gate end

"""
    struct UnitaryGate{F <: Function} <: Gate

Represents a unitary quantum gate.

This struct is parametric, allowing it to hold a type-stable function that generates the corresponding unitary operation of the unitary quantum gate.

# Type parameters
- `F <: Function`: The concrete type of the factory function that creates the unitary operation.

# Fields
- `unitaryOperationFactory::F` — A callable object (typically a closure that captures all the necessary information) that, when invoked, returns a `UnitaryOperation`.
- `qubits::Qubits` — The qubits this gate operates on.
- `name::String` — The name of the gate.

# See also
- [`Gate`](@ref)
- [`UnitaryOperation`](@ref)
- [`Qubits`](@ref)

# Example
Add a Hadamard gate to a 3-qubit quantum circuit at qubit position 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> hGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")
```
"""
struct UnitaryGate{F <: Function} <: Gate
	unitaryOperationFactory::F
	qubits::Qubits
	name::String
end

"""
    qubits(
		gate::Gate
	) -> Qubits

Extracts the list of qubits that a gate is connected to.

# Supported types
- `UnitaryGate`
- `MeasureGate`
- `QuantumChannelGate`

# Returns
- `Qubits` — A vector of qubit indices.

# See also
- [`Gate`](@ref)
- [`Qubits`](@ref)
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)

# Example
Extract the qubits connection of a cnot, measurement and quantum channel gate:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> cnotGate!(qc, 2, 3)
julia> measureGate!(qc, [1, 3])
julia> qubits(getGate(qc, getStep(qc, 1), 1))
2-element Vector{Int64}:
 2
 3
julia> qubits(getGate(qc, getStep(qc, 2), 1))
2-element Vector{Int64}:
 1
 3
```
"""
function qubits(gate::UnitaryGate)
	gate.qubits
end

"""
    Sigmas = Vector{ComplexMatrix}

Alias for a vector of **measurement direction operators**, with one operator per measured qubit.

A **measurement direction operator** is defined on an arbitrary axis on the Bloch sphere and is used to specify the measurement basis. Each **measurement direction operator** can be parameterized by a 3D unit vector ``n = [n_x, n_y, n_z]`` on the Bloch sphere:
```math
    M(n) = n ⋅ σ = n_x ⋅ σ_x + n_y ⋅ σ_y + n_z ⋅ σ_z
```

with ``σ = [σ_x, σ_y, σ_z]`` the vector of Pauli matrices associated with the X, Y, and Z axes, respectively.

Alternatively, if the unit vector ``n`` is expressed in terms of spherical polar angles ``θ`` and ``φ``, the operator takes the form:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

Both forms define the same 2×2 Hermitian operator, specifying a measurement axis on the Bloch sphere as a linear combination of Pauli matrices ``σ_x``, ``σ_y`` and ``σ_z``.

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaY`](@ref)
- [`sigmaZ`](@ref)
- [`sigmaN`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

"""
Sigmas = Vector{ComplexMatrix}

"""
    struct KrausOperator

Represents a single Kraus operator in matrix form, used in quantum channels.

# Fields
- `E::Matrix{ComplexF64}` — Matrix representation of the Kraus operator ``E``.
- `label::String` — A label describing the operator (e.g., measurement outcome).
"""
struct KrausOperator
	E::ComplexMatrix
	label::String
end

"""
	KrausOperators = Vector{KrausOperator}

Alias for a vector of `KrausOperator`s.

# See also
- [`KrausOperator`](@ref)

# Example
Generate the Kraus operators for a projective measurement along direction ``σ_z`` on qubit 1 in a 2-qubit system:
```julia-repl
julia> krausOperators = generateKrausOperatorsForPVMMeasurement(2, [1], [sigmaZ()])
julia> krausOperators[1].E
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
julia> krausOperators[1].label
"0*"
julia> krausOperators[2].E
4×4 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
julia> krausOperators[2].label
"1*"
```
"""
KrausOperators = Vector{KrausOperator}

"""
    @enum MeasureGateType PVM POVM

Supported measurement gate types:

- `PVM`: **Projective Valued Measurement**  
  Represents a standard quantum measurement in which each measurement outcome corresponds to a projection operator. These measurements are described by a set of orthogonal projectors that sum to the identity.

- `POVM`: **Positive Operator-Valued Measure**
  Represents a generalized measurement that can include non-orthogonal measurement elements. These are described by a set of positive semi-definite operators that also sum to the identity, allowing for more general forms of measurement beyond projective ones.

# See also
- [`measureGate!`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
"""
@enum MeasureGateType PVM POVM

"""
    struct MeasureGate <: Gate

Represents a quantum measurement gate.

# Fields
- `measureGateType::MeasureGateType` - Measurenent gate type (PVM or POVM).
- `qubits::Qubits` — Qubits being measured.
- `sigmas::Sigmas` — Measurement direction operators, one for each qubit being measured. Used for **projective measurement (PVM)**.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the **generalized measurement (POVM)**.
- `forgetOutcome::Bool` — If `true`, the measurement result is not revealed (equivalent to taking the partial trace). If `false`, it is observed.
  - `false` - Simulates measurement with revealing outcome by collapsing the state to one of the possible outcomes.
  - `true` - Simulates measurement without revealing outcome by creating a mixed state over all possible outcomes.

# Constructors
- `MeasureGate(qubits::Qubits, sigmas::Sigmas, forgetOutcome::Bool)` - **Projective measurement (PVM)**.
- `MeasureGate(qubits::Qubits, krausOperators::KrausOperators, forgetOutcome::Bool)` - **Generalized measurement (POVM)**.

# See also
- [`Gate`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Qubits`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)

# Examples

**Example 1: Add a PVM measurement gate to a single-qubit quantum circuit on qubit 1 with measurement direction operator ``σ_z`` where the outcome is revealed to the observer**
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> measureGate!(qc, [1], [sigmaZ()], forgetOutcome = false)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1], Matrix{ComplexF64}[[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 0.0im]], nothing, false)
```

**Example 2: Add a POVM measurement gate to a single-qubit quantum circuit using Kraus operators for a projective measurement along direction ``σ_z`` where the outcome is revealed to the observer**
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> measureGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]), forgetOutcome = false)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.POVM, [1], nothing, KrausOperator[KrausOperator(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], "0"), KrausOperator(ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], "1")], false)
```
"""
struct MeasureGate <: Gate
	measureGateType::MeasureGateType
	qubits::Qubits
	sigmas::Union{Nothing, Sigmas}
	krausOperators::Union{Nothing, KrausOperators}
	forgetOutcome::Bool
    function MeasureGate(qubits::Qubits, sigmas::Sigmas, forgetOutcome::Bool)
        new(PVM, qubits, sigmas, nothing, forgetOutcome)
    end
    function MeasureGate(qubits::Qubits, krausOperators::KrausOperators, forgetOutcome::Bool)
        new(POVM, qubits, nothing, krausOperators, forgetOutcome)
    end
end

function qubits(gate::MeasureGate)
	gate.qubits
end

"""
    struct QuantumChannelGate <: Gate

Represents a **quantum channel gate**, also known as a *completely positive trace-preserving (CPTP) map*.

A quantum channel models a **superoperator** — a general type of quantum transformation that includes noise, irreversible processes, decoherence, collapse due to (partial) measurement, or general open system dynamics.

# Fields
- `qubits::Qubits` — Qubits indices on which the quantum channel acts.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the quantum channel.

# See also
- [`Gate`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`Qubits`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)

# Example
Add a **quantum channel gate** to a single-qubit circuit using Kraus operators for a projective measurement along direction ``σ_z`` with the outcome not revealed:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
QuantumChannelGate([1], KrausOperator[KrausOperator(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], "0"), KrausOperator(ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], "1")])
```
"""
struct QuantumChannelGate <: Gate
	qubits::Qubits
	krausOperators::KrausOperators
end

function qubits(gate::QuantumChannelGate)
	gate.qubits
end


"""
    struct BernoulliGate{F <: Function} <: Gate

Represents a Bernoulli quantum gate.

This struct is parametric, allowing it to hold a type-stable function that generates the corresponding unitary operation of the Bernoulli quantum gate.

# Type parameters
- `F <: Function`: The concrete type of the factory function that creates the unitary operation.

# Fields
- `unitaryOperationFactory::F` — A callable object (typically a closure that captures all the necessary information) that, when invoked, returns a `UnitaryOperation`.
- `qubits::Qubits` — The qubits this gate operates on.
- `probability::Float64` — The probability of applying the unitary operation.

# See also
- [`Gate`](@ref)
- [`UnitaryOperation`](@ref)
- [`BernoulliOperation`](@ref)
- [`Qubits`](@ref)

# Example
Add a Bernoulli gate for a Pauli-X operation with probability 0.3 to a 3-qubit quantum circuit at qubit position 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> bernoulliGate!(qc, [2], 0.3, createSingleQubitOperationX())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
BernoulliGate{QubiSim.var"#129#130"{UnitaryOperation}}(QubiSim.var"#129#130"{UnitaryOperation}(UnitaryOperation(ComplexF64[0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im])), [2], 0.3)
```
"""
struct BernoulliGate{F <: Function} <: Gate
	unitaryOperationFactory::F
	qubits::Qubits
	probability::Float64
end

function qubits(gate::BernoulliGate)
	gate.qubits
end

"""
    struct StochasticGate{F <: Function} <: Gate

Represents a stochastic quantum gate.

This struct is parametric, allowing it to hold a type-stable function that generates the corresponding unitary operation of the stochastic quantum gate.

# Type parameters
- `F <: Function`: The concrete type of the factory function that creates the unitary operation.

# Fields
- `unitaryOperationFactory::F` — A callable object (typically a closure that captures all the necessary information) that, when invoked, returns a `UnitaryOperation`.
- `qubits::Qubits` — The qubits this gate operates on.
- `accuracyCheckForUnitarity::Float64` — Numerical tolerance for unitarity validation.

# See also
- [`Gate`](@ref)
- [`UnitaryOperation`](@ref)
- [`StochasticOperation`](@ref)
- [`Qubits`](@ref)

# Example
Add a stochastic gate for a U₁(λ) unitary operation with random phase rotation angle ``λ`` to a 3-qubit quantum circuit at qubit position 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> stochasticGate!(qc, [2], () -> createSingleQubitOperationU1(rand()))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
StochasticGate{var"#11#12"}(var"#11#12"(), [2], 2.220446049250313e-15)
```
"""
struct StochasticGate{F <: Function} <: Gate
	unitaryOperationFactory::F
	qubits::Qubits
	accuracyCheckForUnitarity::Float64
end

function qubits(gate::StochasticGate)
	gate.qubits
end

"""
	abstract type Step

Abstract type representing a generic step in a quantum circuit.
"""
abstract type Step end

"""
    struct UnitaryStep <: Step

Represents a unitary step in a quantum circuit, composed of a vector of unitary gates.

# Fields
- `gates::Vector{UnitaryGate}` — A list of unitary gates applied in this `Step`.

# See also
- [`Step`](@ref)
- [`UnitaryGate`](@ref)
"""
struct UnitaryStep <: Step
    gates::Vector{UnitaryGate}
end

"""
    struct MeasurementStep <: Step

Represents a measurement step in a quantum circuit, consisting of a single measurement gate.

# Fields
- `gate::MeasureGate` — The measurement gate used in this `Step`.

# See also
- [`Step`](@ref)
- [`MeasureGate`](@ref)
"""
struct MeasurementStep <: Step
    gate::MeasureGate
end

"""
    struct QuantumChannelStep <: Step

Represents a quantum channel step in a quantum circuit, consisting of a single quantum channel gate.

# Fields
- `gate::QuantumChannelGate` — The quantum channel gate used in this `Step`.

# See also
- [`Step`](@ref)
- [`QuantumChannelGate`](@ref)
"""
struct QuantumChannelStep <: Step
    gate::QuantumChannelGate
end

"""
    struct BernoulliStep <: Step

Represents a Bernoulli step in a quantum circuit, consisting of a single Bernoulli gate.

# Fields
- `gate::BernoulliGate` — The Bernoulli gate used in this `Step`.

# See also
- [`Step`](@ref)
- [`BernoulliGate`](@ref)
"""
struct BernoulliStep <: Step
    gate::BernoulliGate
end

"""
    struct StochasticStep <: Step

Represents a stochastic step in a quantum circuit, consisting of a single stochastic gate.

# Fields
- `gate::StochasticGate` — The stochastic gate used in this `Step`.

# See also
- [`Step`](@ref)
- [`StochasticGate`](@ref)
"""
struct StochasticStep <: Step
    gate::StochasticGate
end

"""
	Circuit = Vector{Step}

Type alias for a quantum circuit, defined as a vector of `Step` instances.

# See also
- [`Step`](@ref)
"""
Circuit = Vector{Step}

"""
    struct QuantumCircuit{Z, I}

Represents a full quantum circuit, including unitary gates, measurements, and quantum channels.

This struct is **compile-time specialized** on the type parameters `Z` and `I` (that encode the qubit indexing settings), ensuring optimized and branch-free code generation.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Fields
- `settings::Settings{Z, I}` — Qubit indexing settings.
- `numberOfQubits::Int` — The total number of qubits in the circuit.
- `numberOfStepsOnQubits::Vector{Int}` — Tracks the number of steps applied to each qubit.
- `circuit::Circuit` — A sequence of circuit steps, each being a subtype of `Step` (`UnitaryStep`, `MeasurementStep`, `QuantumChannelStep`).

# See also
- [`Circuit`](@ref)
- [`Step`](@ref)
- [`Settings`](@ref)
- [`UnitaryStep`](@ref)
- [`MeasurementStep`](@ref)
- [`QuantumChannelStep`](@ref)
"""
mutable struct QuantumCircuit{Z, I}
	settings::Settings{Z, I}
	numberOfQubits::Int
	numberOfStepsOnQubits::Vector{Int}
	circuit::Circuit
end

"""
    getStep(
		quantumCircuit::QuantumCircuit, 
		stepId::Int
	) -> Step

Extracts a specific step from a quantum circuit by its ID.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit from which to extract the step.
- `stepId::Int` — The step index (1-based or 0-based depending on `Settings`).

# Returns
- A `Step` from the circuit, such as a `UnitaryStep`, `MeasurementStep`, or `QuantumChannelStep`.

# See also
- [`QuantumCircuit`](@ref)
- [`Step`](@ref)
- [`Settings`](@ref)
- [`UnitaryStep`](@ref)
- [`MeasurementStep`](@ref)
- [`QuantumChannelStep`](@ref)

# Example
Extract the first and second steps from a quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> hGate!(qc, 2)
julia> cnotGate!(qc, 1, 2)
julia> getStep(qc, 1)
UnitaryStep(UnitaryGate[
	UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H"), 
	UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")
])
julia> getStep(qc, 2)
UnitaryStep(UnitaryGate[
	UnitaryGate{QubiSim.var"#41#42"}(QubiSim.var"#41#42"(), [1, 2], "CNOT")
])
```
"""
function getStep(quantumCircuit::QuantumCircuit, stepId::Int)
	return quantumCircuit.circuit[convertToOneBasedNumbering(quantumCircuit.settings, stepId)]
end

"""
    getGate(
		quantumCircuit::QuantumCircuit, 
		step, 
		gateId::Int
	) -> Gate

Retrieves a specific gate from a step in the quantum circuit.

# Arguments
- `quantumCircuit::QuantumCircuit` — The parent circuit.
- `step` — A circuit step (e.g. `UnitaryStep`, `MeasurementStep`, `QuantumChannelStep`, `BernoulliStep`, `StochasticStep`).
- `gateId::Int` — The index of the gate within the step (1-based or 0-based depending on `Settings`).

# Returns
- The specified `Gate` (e.g. `UnitaryGate`, `MeasureGate`, `QuantumChannelGate`, `BernoulliGate`, or `StochasticGate`).

# Notes
- `UnitaryStep`s can contain multiple gates.
- `MeasurementStep`, `QuantumChannelStep`, `BernoulliStep` and `StochasticStep` each contain only one gate; `gateId` must be 1.

# See also
- [`QuantumCircuit`](@ref)
- [`Gate`](@ref)
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)
- [`BernoulliGate`](@ref)
- [`StochasticGate`](@ref)
- [`Settings`](@ref)
- [`UnitaryStep`](@ref)
- [`MeasurementStep`](@ref)
- [`QuantumChannelStep`](@ref)
- [`BernoulliStep`](@ref)
- [`StochasticStep`](@ref)

# Example
Extract the gates from a quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1)
UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")
julia> getGate(qc, getStep(qc, 1), 2)
UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")
julia> getGate(qc, getStep(qc, 2), 1)
UnitaryGate{QubiSim.var"#41#42"}(QubiSim.var"#41#42"(), [1, 2], "CNOT")
```
"""
function getGate(quantumCircuit::QuantumCircuit, step::UnitaryStep, gateId::Int)
	return step.gates[convertToOneBasedNumbering(quantumCircuit.settings, gateId)]
end
function getGate(quantumCircuit::QuantumCircuit, step::MeasurementStep, gateId::Int)
	if convertToOneBasedNumbering(quantumCircuit.settings, gateId) == 1
		return step.gate
	else
		throw("The gateId=$(gateId) and cannot be greater than 1 since a MeasurementStep contains only 1 gate.")
	end
end
function getGate(quantumCircuit::QuantumCircuit, step::QuantumChannelStep, gateId::Int)
	if convertToOneBasedNumbering(quantumCircuit.settings, gateId) == 1
		return step.gate
	else
		throw("The gateId=$(gateId) and cannot be greater than 1 since a QuantumChannelStep contains only 1 gate.")
	end
end
function getGate(quantumCircuit::QuantumCircuit, step::BernoulliStep, gateId::Int)
	if convertToOneBasedNumbering(quantumCircuit.settings, gateId) == 1
		return step.gate
	else
		throw("The gateId=$(gateId) and cannot be greater than 1 since a BernoulliStep contains only 1 gate.")
	end
end
function getGate(quantumCircuit::QuantumCircuit, step::StochasticStep, gateId::Int)
	if convertToOneBasedNumbering(quantumCircuit.settings, gateId) == 1
		return step.gate
	else
		throw("The gateId=$(gateId) and cannot be greater than 1 since a StochasticStep contains only 1 gate.")
	end
end
