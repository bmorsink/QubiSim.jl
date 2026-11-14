############################################
# QuantumProgramDefinitions

"""
	abstract type Operation

Abstract base type for all quantum operations
"""
abstract type Operation end

"""
	Program = Vector{Operation}

Type alias for a quantum program, defined as a sequence (vector) of quantum `Operation` instances.

# See also
- [`Operation`](@ref)
"""
Program = Vector{Operation}

""" 
    struct QuantumProgram{Z, I}

Represents a compiled quantum circuit as a sequence of operations.

This struct is **compile-time specialized** on the type parameters `Z` and `I` (that encode the qubit indexing settings), ensuring optimized and branch-free code generation.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Fields
- `settings::Settings{Z, I}` — Qubit indexing settings.
- `numberOfQubits::Int` — Total number of qubits used in the program.
- `program::Program` — Sequence of `Operation`s applied to the qubits.

# Subtypes of `Operation` are `UnitaryOperation`, `MeasureOperation`, `MeasureAndForgetOperation`, or `QuantumChannelOperation`.

# See also
- [`Settings`](@ref)
- [`Program`](@ref)
- [`Operation`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)
"""
struct QuantumProgram{Z, I}
	settings::Settings{Z, I}
	numberOfQubits::Int
	program::Program
end

""" 
    struct UnitaryOperation <: Operation

A unitary operation represented by a unitary matrix.

# Fields
- `U::ComplexMatrix` — Unitary matrix ``U`` defining the operation.

# Quantum state evolution
- **Quantum vector state**:  
    - A pure state ``|ψ⟩`` evolves as: ``|ψ⟩ → U·|ψ⟩``
- **Quantum density state**:
    - A density operator ``ρ`` evolves as: ``ρ → U·ρ·U^†``

# See also
- [`UnitaryGate`](@ref)
- [`Operation`](@ref)

# Example
Extract the unitary operation (its matrix) of a Hadamard gate added to a single-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> qp = compileQuantumCircuit(qc)
julia> getOperation(qp, 1).U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
struct UnitaryOperation <: Operation
	U::ComplexMatrix
end

UnitaryOperation(U::AbstractMatrix{<:Number}) = UnitaryOperation(ComplexF64.(U))

""" 
    struct MeasureOperation <: Operation

Represents a quantum measurement operation on a subset of qubits. The measurement result is revealed to the observer.
Creates a post-measurement pure state by collapsing the state to one of the possible outcomes.

# Fields
- `krausOperators::KrausOperators` — A collection of Kraus operators ``{Eₖ}`` that defines the measurement process.
- `measurementOperators::Vector{ComplexMatrix}` — A collection of measurement operators ``{Mₖ}`` corresponding to the Kraus operators.
- `labels::Vector{String}` — A vector of binary string labels corresponding to each outcome (e.g., `"00"`, `"01"`, etc.).

# Quantum state evolution under a generalized measurement (POVM)

Let ``{Mₖ}`` be measurement operators satisfying ``∑ₖ Mₖ^†·Mₖ = I``.
Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
- **Quantum vector state**:
    - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·|ψ⟩⟨ψ|)``
    - Post-measurement state for outcome ``k``: ``|ψ⟩ → Mₖ·|ψ⟩/√pₖ``
- **Quantum density state**:
    - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·ρ)``
    - Post-measurement state for outcome ``k``: ``ρ → Eₖ·ρ·Eₖ^†/pₖ``

# Quantum state evolution under a projective measurement (PVM)

A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``.
Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`.
Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

# See also
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Operation`](@ref)

# Example
Extract the measurement operation of a measurement gate added to a 2-qubit quantum circuit at qubit 1 with measurement direction operator ``σ_z`` that reveals the result to the observer:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> measureGate!(qc, [1], [sigmaZ()], forgetOutcome=false)
julia> qp = compileQuantumCircuit(qc)
julia> typeof(getOperation(qp, 1))
MeasureOperation
julia> getOperation(qp, 1).krausOperators[1].E
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
julia> getOperation(qp, 1).krausOperators[1].label
"0*"
julia> getOperation(qp, 1).krausOperators[2].E
4×4 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
julia> getOperation(qp, 1).krausOperators[2].label
"1*"
```
"""
struct MeasureOperation <: Operation
	krausOperators::KrausOperators
	measurementOperators::Vector{ComplexMatrix}
	labels::Vector{String}
end

""" 
    struct MeasureAndForgetOperation <: Operation

Represents a quantum measurement operation on a subset of qubits. The measurement result is NOT revealed to the observer (equivalent to tracing out the subsystem).
Creates a post-measurement mixed state of all possible outcomes.

# Fields
- `krausOperators::KrausOperators` — A collection of Kraus operators ``{Eₖ}`` that defines the measurement process.
- `measurementOperators::Vector{ComplexMatrix}` — A collection of measurement operators ``{Mₖ}`` corresponding to the Kraus operators.
- `labels::Vector{String}` — A vector of binary string labels corresponding to each outcome (e.g., `"00"`, `"01"`, etc.).

# Quantum state evolution under a generalized measurement (POVM)

Let ``{Mₖ}`` be measurement operators satisfying ``∑ₖ Mₖ^†·Mₖ = I``.
Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
- **Quantum vector state**:
    - Combination of a MeasureAndForgetOperation with a quantum vector state is not possible...
- **Quantum density state**:
    - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·ρ)``
    - Post-measurement state (outcome not observed): ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

# Quantum state evolution under a projective measurement (PVM)

A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``.
Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`.
Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

# See also
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`MeasureOperation`](@ref)
- [`Operation`](@ref)

# Example
Extract the measurement operation of a measurement gate added to a 2-qubit quantum circuit at qubit 1 with measurement direction operator ``σ_z`` that does NOT reveal the result to the observer:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> measureGate!(qc, [1], [sigmaZ()], forgetOutcome=true)
julia> qp = compileQuantumCircuit(qc)
julia> typeof(getOperation(qp, 1))
MeasureAndForgetOperation
julia> getOperation(qp, 1).krausOperators[1].E
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
julia> getOperation(qp, 1).krausOperators[1].label
"0*"
julia> getOperation(qp, 1).krausOperators[2].E
4×4 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
julia> getOperation(qp, 1).krausOperators[2].label
"1*"
```
"""
struct MeasureAndForgetOperation <: Operation
	krausOperators::KrausOperators
	measurementOperators::Vector{ComplexMatrix}
	labels::Vector{String}
end

""" 
    struct QuantumChannelOperation <: Operation

Represents a quantum channel — a completely positive, trace-preserving (CPTP) map — applied to one or more qubits.

# Fields
- `krausOperators::KrausOperators` — A collection of Kraus operators ``{Eₖ}`` that defines the quantum channel.

# Quantum state evolution under a quantum channel
- **Quantum density state**:
    - A density operator ``ρ`` evolves as: ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

# See also
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Operation`](@ref)

# Example
Extract the quantum channel operation (its Kraus operators) for an ignorant measurement (using a quantum channel gate) with measurement direction operator ``σ_z`` added to a single-qubit quantum circuit.
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]))
julia> qp = compileQuantumCircuit(qc)
julia> getOperation(qp, 1).krausOperators
2-element Vector{KrausOperator}:
 KrausOperator(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], "0")
 KrausOperator(ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], "1")
```
"""
struct QuantumChannelOperation <: Operation
	krausOperators::KrausOperators
end

""" 
    struct BernoulliOperation <: Operation

A Bernoulli operation represented by a unitary matrix and the probability of applying it.

# Fields
- `U::ComplexMatrix` — Unitary matrix ``U`` defining the operation.
- `probability::Float64` — The probability of applying the operation.

# Quantum state evolution
- **Quantum vector state**:
    - A pure state ``|ψ⟩`` evolves as: ``|ψ⟩ → U·|ψ⟩`` with probability `probability` while with probability `1-probability` the state is left unchanged ``|ψ⟩ → |ψ⟩``.
- **Quantum density state**:
    - A density operator ``ρ`` evolves as: ``ρ → U·ρ·U^†`` with probability `probability` while with probability `1-probability` the state is left unchanged ``ρ → ρ``.

# See also
- [`BernoulliGate`](@ref)
- [`Operation`](@ref)

# Example
Extract the Bernoulli operation of a Bernoulli gate for a Pauli-X operation with probability 0.3 added to a single-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> bernoulliGate!(qc, [1], 0.3, createSingleQubitOperationX())
julia> qp = compileQuantumCircuit(qc)
julia> getOperation(qp, 1)
BernoulliOperation(ComplexF64[0.0 + 0.0im 1.0 - 0.0im; 1.0 + 0.0im 0.0 + 0.0im], 0.3)
```
"""
struct BernoulliOperation <: Operation
	U::ComplexMatrix
	probability::Float64
end

""" 
    struct StochasticOperation <: Operation

A stochastic operation represented by a function that generates the unitary operation.

# Type parameters
- `F <: Function`: The concrete type of the factory function that creates the unitary operation.

# Fields
- `unitaryOperationFactory::F` — A callable object (typically a closure that captures all the necessary information) that, when invoked, returns a `UnitaryOperation`.
- `qubits::Qubits` — The qubits this gate operates on.
- `numberOfQubits::Int` — Total number of qubits used in the program.
- `accuracyCheckForUnitarity::Float64` — Numerical tolerance for unitarity validation.

# Quantum state evolution
During execution of the quantum program, `unitaryOperationFactory` is invoked once per shot to generate a unitary operation ``U``. 
The generated operation is then applied to the quantum state as follows:
- **Quantum vector state**:  
    - A pure state ``|ψ⟩`` evolves as: ``|ψ⟩ → U·|ψ⟩``
- **Quantum density state**:
    - A density operator ``ρ`` evolves as: ``ρ → U·ρ·U^†``

# See also
- [`StochasticGate`](@ref)
- [`Operation`](@ref)

# Example
Extract the stochastic operation of a stochastic gate for a U₁(λ) unitary operation with random phase rotation angle ``λ`` added to a single-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> stochasticGate!(qc, [1], () -> createSingleQubitOperationU1(rand()))
julia> qp = compileQuantumCircuit(qc)
julia> getOperation(qp, 1)
StochasticOperation{var"#13#14"}(var"#13#14"(), [1], 1, 2.220446049250313e-15)
```
"""
struct StochasticOperation{F <: Function} <: Operation
	unitaryOperationFactory::F
	qubits::Qubits
	numberOfQubits::Int
	accuracyCheckForUnitarity::Float64
end

"""
    getOperation(
		quantumProgram::QuantumProgram, 
		operationId::Int
	) -> Operation

Returns the `operationId`-th `Operation` from a `QuantumProgram`.

The `Operation` is retrieved using the program’s internal indexing scheme as defined by its `Settings`.

# Arguments
- `quantumProgram::QuantumProgram` — The compiled quantum circuit containing a sequence of operations.
- `operationId::Int` — The index of the operation to retrieve (1-based or 0-based depending on `settings`).

# Returns
- `Operation` — A subtype of `Operation`, such as a `UnitaryOperation`, `MeasureOperation`, `MeasureAndForgetOperation`, `QuantumChannelOperation`, `BernoulliOperation`, or `StochasticOperation`.

# See also
- [`QuantumProgram`](@ref)
- [`Operation`](@ref)
- [`Settings`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`BernoulliOperation`](@ref)
- [`StochasticOperation`](@ref)

# Example
Extract the first and second operations from a quantum program created from a single-qubit circuit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> qp = compileQuantumCircuit(qc)
julia> getOperation(qp, 1)
UnitaryOperation(ComplexF64[0.7071067811865476 + 0.0im 0.7071067811865475 + 0.0im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im])
julia> getOperation(qp, 2)
UnitaryOperation(ComplexF64[0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im])
```
"""
function getOperation(quantumProgram::QuantumProgram, operationId::Int)
	return quantumProgram.program[convertToOneBasedNumbering(quantumProgram.settings, operationId)]
end