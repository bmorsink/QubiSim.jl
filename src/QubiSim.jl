module QubiSim

export ComplexVector, ComplexMatrix
export IndexType, circuitIndexBigEndian, circuitIndexLittleEndian, byteIndex
export StateSpace, vector, density
export UseFromMeasurement, useProbability, useOutcome

export Settings, createSettings, zeroBasedNumbering, indexType

export Qubits, Gate, Step, Circuit, Sigmas, Operation, Program, State, KrausOperator, KrausOperators
export QuantumCircuit, UnitaryStep, UnitaryGate, qubits, MeasurementStep, MeasureGateType, MeasureGate, QuantumChannelStep, QuantumChannelGate, getStep, getGate
export QuantumProgram, UnitaryOperation, MeasureOperation, MeasureAndForgetOperation, QuantumChannelOperation, getOperation
export QuantumOutput, Measured, VectorState, DensityState, getState, getShot, getStateShot, convertVectorStateToDensityState

export createQuantumCircuit, createQuantumCircuitWithSettings, concatenateQuantumCircuits, inverseQuantumCircuit
export compileQuantumCircuit, compileToSingleGate

export createByteIndexVector, createInitialQubitState, createDoubleQubitWernerDensityState, createSingleQubitBlochDensityState
export runQuantumProgram, calculateEntropyAndPurity
export calculateCorrelation, calculateExpectationValueOfProductOfMeasuredQubits

export createToggleSwapList, tensorProduct, partialTrace, fidelity

export barrier!
export unitaryUGate!
export measureGate!, sigmaX, sigmaY, sigmaZ, sigmaN, pauliX, pauliY, pauliZ, generateKrausOperatorsForPVMMeasurement
export generateKrausOperatorsForDepolarizingChannel, generateKrausOperatorsForPhaseDampingChannel, generateKrausOperatorsForAmplitudeDampingChannel
export quantumChannelGate!
export u1Gate!, createSingleQubitOperationU1
export u2Gate!, createSingleQubitOperationU2
export u3Gate!, createSingleQubitOperationU3
export xGate!, createSingleQubitOperationX
export yGate!, createSingleQubitOperationY
export zGate!, createSingleQubitOperationZ
export hGate!, createSingleQubitOperationH
export idGate!, createSingleQubitOperationId
export rxGate!, createSingleQubitOperationRx
export ryGate!, createSingleQubitOperationRy
export rzGate!, createSingleQubitOperationRz
export rotationGate!, createSingleQubitOperationRotation
export projectionGate!, createNQubitOperationProjection
export reflectionGate!, createNQubitOperationReflection
export tGate!, createSingleQubitOperationT
export tdGate!, createSingleQubitOperationTd
export sGate!, createSingleQubitOperationS
export sdGate!, createSingleQubitOperationSd
export cnotGate!, createDoubleQubitOperationCNOT
export cnotReverseGate!, createDoubleQubitOperationCNOTReverse
export swapGate!, createDoubleQubitOperationSWAP
export phaseGate!, createDoubleQubitOperationPhase
export expHGate!, createNQubitOperationExpH
export qftGate!, createNQubitOperationQFT, createNQubitQFTQuantumCircuit
export iqftGate!, createNQubitOperationIQFT, createNQubitIQFTQuantumCircuit
export controlledUGate!, createNQubitOperationControlledU
export controlledHGate!, createDoubleQubitOperationControlledH
export controlledXGate!, createDoubleQubitOperationControlledX
export controlledYGate!, createDoubleQubitOperationControlledY
export controlledZGate!, createDoubleQubitOperationControlledZ
export controlledTGate!, createDoubleQubitOperationControlledT
export controlledTdGate!, createDoubleQubitOperationControlledTd
export controlledSGate!, createDoubleQubitOperationControlledS
export controlledSdGate!, createDoubleQubitOperationControlledSd
export qpeGate!, createNQubitOperationQPE, createNQubitQPEQuantumCircuit

export probeMeasureOutcome, probeMeasureProbability, probeStateProbability, probeStateMultiBlochVector

using LinearAlgebra, PlotlyJS, StatsBase
import Base: broadcastable

######################
# QuantumConfiguration

""" 
	@enum IndexType circuitIndexBigEndian circuitIndexLittleEndian byteIndex

Enumeration defining the qubit indexing strategy used in a quantum circuit.

# Variants

- `circuitIndexBigEndian`: 
  - Qubits are indexed top-to-bottom in the circuit view (1 to n).
  - This corresponds to left-to-right (1 to n) in the byte representation.

- `circuitIndexLittleEndian`: 
  - Qubits are indexed top-to-bottom (1 to n) in the circuit view.
  - This corresponds to right-to-left (n to 1) in the byte representation.

- `byteIndex`: 
  - Qubits are indexed left-to-right (1 to n) directly in byte form.

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# See also
- [`Settings`](@ref)
"""
@enum IndexType circuitIndexBigEndian circuitIndexLittleEndian byteIndex

""" 
	struct Settings{Z, I}

Compile-time structure defining **global indexing settings** for qubits in a quantum circuit.  
This type encodes its parameters as type-level values for efficient specialization and dispatch.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Constructors
- `Settings(zeroBasedNumbering::Bool, indexType::IndexType)` - Creates a specialized Settings{Z, I} instance, where both parameters `Z` and `I` are encoded at compile time.

# Accessors
- `zeroBasedNumbering(settings::Settings) -> Bool` - Returns whether qubit indices start from 0 or 1:
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1.

- `indexType(settings::Settings) -> IndexType` - Returns the index convention used:
  - Determines how qubits are mapped between circuit view and byte representation. 
  - Can be `circuitIndexBigEndian`, `circuitIndexLittleEndian`, or `byteIndex`.

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# See also
- [`IndexType`](@ref)

# Example
```julia-repl
julia> s = Settings(true, circuitIndexBigEndian)
Settings{true, circuitIndexBigEndian}()
julia> zeroBasedNumbering(s)
true
julia> indexType(s)
circuitIndexBigEndian::IndexType = 0
```
"""
struct Settings{Z, I} end
Settings(zeroBasedNumbering::Bool, indexType::IndexType) = Settings{zeroBasedNumbering, indexType}()
zeroBasedNumbering(::Settings{Z, I}) where {Z, I} = Z
indexType(::Settings{Z, I}) where {Z, I} = I

# """
#     convertToOneBasedNumbering(
# 		settings::Settings, 
# 		index
# 	) -> Int

# Converts an index to **1-based numbering** based on the provided settings.

# This function is **compile-time specialized** on the `zeroBasedNumbering` parameter of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `zeroBasedNumbering = true`, the function **adds 1** to the given `index`.
# - If `zeroBasedNumbering = false`, the `index` is **returned unchanged**.

# # See also
# - [`Settings`](@ref)
# """
convertToOneBasedNumbering(::Settings{true, I}, index) where {I} = index.+1
convertToOneBasedNumbering(::Settings{false, I}, index) where {I} = index

# """
#     convertFromOneToSettingsBasedNumbering(
# 		settings::Settings, 
# 		index
# 	) -> Int

# Converts a **1-based index** to the numbering convention specified in the given `settings`.

# This function is **compile-time specialized** on the `zeroBasedNumbering` parameter of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `zeroBasedNumbering = true`, the function **subtracts 1** from the given `index`.
# - If `zeroBasedNumbering = false`, the `index` is **returned unchanged**.

# # See also
# - [`Settings`](@ref)
# """
convertFromOneToSettingsBasedNumbering(::Settings{true, I}, index) where {I} = index.-1
convertFromOneToSettingsBasedNumbering(::Settings{false, I}, index) where {I} = index

# """
#     convertToByteIndex(
# 		settings::Settings, 
# 		numberOfQubits::Int, 
# 		index
# 	) -> Int

# Converts a **qubit index** to its corresponding **byte index**, taking into account both the **endianness** and the **indexing base** defined in the given `settings`.

# This function is **compile-time specialized** on the `indexType` and `zeroBasedNumbering` parameters of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `indexType = circuitIndexBigEndian` or `indexType = byteIndex`, the function returns the same (but one-based converted) `index`.
# - If `indexType = circuitIndexLittleEndian`, the **flipped** one-based converted `index` relative to the `numberOfQubits` is returned.

# # See also
# - [`Settings`](@ref)
# ```
# """
convertToByteIndex(settings::Settings{Z, I}, numberOfQubits::Int, index) where {Z, I} = convertToOneBasedNumbering(settings, index)
convertToByteIndex(settings::Settings{Z, circuitIndexLittleEndian}, numberOfQubits::Int, index) where {Z} = numberOfQubits.-convertToOneBasedNumbering(settings, index).+1

"""
	ComplexVector = Vector{ComplexF64}

Alias for a vector of complex values.
"""
ComplexVector = Vector{ComplexF64}

"""
	ComplexMatrix = Matrix{ComplexF64}

Alias for a matrix of complex values.
"""
ComplexMatrix = Matrix{ComplexF64}

######################
# QuantumCircuit

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
- `step` — A circuit step (e.g. `UnitaryStep`, `MeasurementStep`, `QuantumChannelStep`).
- `gateId::Int` — The index of the gate within the step (1-based or 0-based depending on `Settings`).

# Returns
- The specified `Gate` (e.g. `UnitaryGate`, `MeasureGate`, or `QuantumChannelGate`).

# Notes
- `UnitaryStep`s can contain multiple gates.
- `MeasurementStep` and `QuantumChannelStep` each contain only one gate; `gateId` must be 1.

# See also
- [`QuantumCircuit`](@ref)
- [`Gate`](@ref)
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)
- [`Settings`](@ref)
- [`UnitaryStep`](@ref)
- [`MeasurementStep`](@ref)
- [`QuantumChannelStep`](@ref)

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

######################
# QuantumProgram

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
- `Operation` — A subtype of `Operation`, such as a `UnitaryOperation`, `MeasureOperation`, `MeasureAndForgetOperation`, or `QuantumChannelOperation`.

# See also
- [`QuantumProgram`](@ref)
- [`Operation`](@ref)
- [`Settings`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)

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

######################
# QuantumOutput

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

######################

function addGate!(quantumCircuit::QuantumCircuit, gate::UnitaryGate)
	maximumNumberOfStepsOnQubits = maximum(quantumCircuit.numberOfStepsOnQubits[qubits(gate)])

	if length(quantumCircuit.circuit) < maximumNumberOfStepsOnQubits
		step = UnitaryStep([gate])
		push!(quantumCircuit.circuit, step)
	else
		push!(quantumCircuit.circuit[maximumNumberOfStepsOnQubits].gates, gate)
	end
	
	for qubit in qubits(gate)
		quantumCircuit.numberOfStepsOnQubits[qubit] = maximumNumberOfStepsOnQubits + 1
	end
end
function addGate!(quantumCircuit::QuantumCircuit, gate::QuantumChannelGate)
	push!(quantumCircuit.circuit, QuantumChannelStep(gate))
	quantumCircuit.numberOfStepsOnQubits = (length(quantumCircuit.circuit)+1)*ones(Int, quantumCircuit.numberOfQubits)
end
function addGate!(quantumCircuit::QuantumCircuit, gate::MeasureGate)
	push!(quantumCircuit.circuit, MeasurementStep(gate))
	quantumCircuit.numberOfStepsOnQubits = (length(quantumCircuit.circuit)+1)*ones(Int, quantumCircuit.numberOfQubits)
end

""" 
    barrier!(
		quantumCircuit::QuantumCircuit
	) -> Nothing

Inserts a **barrier** into the `quantumCircuit`. A barrier acts as a synchronization point, ensuring that no gate operations are reordered across it during circuit optimization or scheduling.

# Arguments
- `quantumCircuit::QuantumCircuit` - The quantum circuit to which the barrier will be added.

# Description
A barrier prevents transformations or gate optimizations from reordering operations across it, effectively separating different stages or logical steps in a quantum circuit.

This is particularly useful for:
- Visual clarity when visualizing quantum circuits
- Debugging or ensuring measurement occurs after specific gates
- Segmenting logical layers of gates

Internally, it updates the `numberOfStepsOnQubits` vector so that all qubits align to the maximum step depth.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the barrier in-place.

# See also
- [`QuantumCircuit`](@ref)

# Example
Add a barrier to a 2-qubit quantum circuit:
```julia-repl
julia> qc=createQuantumCircuit(2)
julia> barrier!(qc)
```
"""
function barrier!(quantumCircuit::QuantumCircuit)
	maximumNumberOfStepsOnQubits = maximum(quantumCircuit.numberOfStepsOnQubits)
	quantumCircuit.numberOfStepsOnQubits = maximumNumberOfStepsOnQubits * ones(Int64, quantumCircuit.numberOfQubits)
	nothing
end

function createSingleQubitOperationU(theta::Float64, phi::Float64, lambda::Float64)
	UnitaryOperation([cos(theta/2) -exp(im*lambda)*sin(theta/2);
             exp(im*phi)*sin(theta/2) exp(im*lambda+im*phi)*cos(theta/2)])
end

function createNQubitOperationId(numberOfQubits::Int)
	UnitaryOperation(Matrix{ComplexF64}(I, 2^numberOfQubits, 2^numberOfQubits))
end

function createNQubitOperatorId(numberOfQubits::Int)
	Matrix{ComplexF64}(I, 2^numberOfQubits, 2^numberOfQubits)
end

#

""" 
	xGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-X gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-X gate** (also known as a quantum NOT gate) performs a unitary operation that flips the state of the qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-X unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationX` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationX`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-X gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc=createQuantumCircuit(3)
julia> xGate!(qc,2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [2], "X")
```
"""
function xGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationX(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"X"
    ))
end

""" 
    createSingleQubitOperationX() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-X unitary operation**.

The **Pauli-X unitary operation** flips the state of the qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

This unitary operation is equivalent to applying the ``U_3(π, 0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-X unitary operation**.

# See also
- [`xGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-X unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationX()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
function createSingleQubitOperationX()
    createSingleQubitOperationU3(pi/1., 0., pi/1.)
end

#

""" 
	yGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-Y gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-Y gate** performs a unitary operation that introduces both a bit and phase flip on a single qubit.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-Y unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationY` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationY`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-Y gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> yGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#7#8"}(QubiSim.var"#7#8"(), [2], "Y")
```
"""
function yGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationY(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Y"
    ))
end

""" 
    createSingleQubitOperationY() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-Y unitary operation**.

The **Pauli-Y unitary operation** flips both the bit and phase of a single qubit.

This unitary operation is equivalent to applying the ``U_3(π, π/2, π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-Y unitary operation**.

# See also
- [`yGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-Y unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationY()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
function createSingleQubitOperationY()
    createSingleQubitOperationU3(pi/1., pi/2., pi/2.)
end

#

""" 
	zGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-Z gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-Z gate** performs a unitary operation that flips the phase of a single qubit.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-Z unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationZ` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationZ`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-Z gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> zGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#9#10"}(QubiSim.var"#9#10"(), [2], "Z")
```
"""
function zGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationZ(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Z"
    ))
end

""" 
    createSingleQubitOperationZ() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-Z unitary operation**.

The **Pauli-Z unitary operation** flips the phase of a single qubit.

This unitary operation is equivalent to applying the ``U_1(π)`` or ``U_3(0, 0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-Z unitary operation**.

# See also
- [`zGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-Z unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationZ()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function createSingleQubitOperationZ()
    createSingleQubitOperationU1(pi/1.)
end

#

""" 
	hGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Hadamard gate** to the `quantumCircuit` on the specified `qubit`.

The **Hadamard gate** performs a unitary operation that creates a superposition of the ``|0⟩`` and ``|1⟩`` states. It performs both a **bit and phase rotation**, mapping computational basis states as follows:

- ``|0⟩ → (|0⟩ + |1⟩) / √2``
- ``|1⟩ → (|0⟩ - |1⟩) / √2``

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Hadamard unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Hadamard gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> hGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")
```
"""
function hGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationH(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"H"
    ))
end

""" 
    createSingleQubitOperationH() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Hadamard unitary operation**.

The **Hadamard unitary operation** creates a superposition of the ``|0⟩`` and ``|1⟩`` states. It performs both a **bit and phase rotation**, mapping computational basis states as follows:

- ``|0⟩ → (|0⟩ + |1⟩) / √2``
- ``|1⟩ → (|0⟩ - |1⟩) / √2``

This unitary operation is equivalent to applying the ``U_2(0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Hadamard unitary operation**.

# See also
- [`hGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Hadamard unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationH()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.707107+0.0im   0.707107+0.0im
0.707107+0.0im  -0.707107+0.0im
```
"""
function createSingleQubitOperationH()
    createSingleQubitOperationU2(0., pi/1.)
end

#

""" 
	idGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **identity gate** to the `quantumCircuit` on the specified `qubit`.

The **identity gate** performs a unitary operation that leaves the qubit state unchanged.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **identity unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationId` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationId`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **identity gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> idGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#13#14"}(QubiSim.var"#13#14"(), [2], "Id")
```
"""
function idGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationId(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Id"
    ))
end

""" 
    createSingleQubitOperationId() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **identity unitary operation**.

The **identity unitary operation** leaves the qubit state unchanged.

This unitary operation is equivalent to applying the ``U_1(0)`` or ``U_3(0, 0, 0)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **identity unitary operation**.

# See also
- [`idGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **identity unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationId()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  1.0+0.0im
```
"""
function createSingleQubitOperationId()
    createSingleQubitOperationU1(0.)
end

#

""" 
	rxGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RX rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RX rotation gate** performs a unitary operation that rotates the qubit state around the **X-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RX rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRx` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRx`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RX rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rxGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#15#16"{Float64}}(QubiSim.var"#15#16"{Float64}(0.5), [2], "RX")
```
"""
function rxGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRx(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RX"
    ))
end

""" 
    createSingleQubitOperationRx(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RX rotation unitary operation**.

The **RX rotation unitary operation** rotates the qubit state around the **X-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RX rotation unitary operation**.

# See also
- [`rxGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RX rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRx(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912-0.0im  0.0-0.247404im
0.0-0.247404im  0.968912+0.0im
```
"""
function createSingleQubitOperationRx(theta::Float64)
    # createSingleQubitOperationU3(theta, -pi/2., pi/2.)
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(pi/2, 0.)))
end

#

""" 
	ryGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RY rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RY rotation gate** performs a unitary operation that rotates the qubit state around the **Y-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RY rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRy` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRy`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RY rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> ryGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#17#18"{Float64}}(QubiSim.var"#17#18"{Float64}(0.5), [2], "RY")
```
"""
function ryGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRy(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RY"
    ))
end

""" 
    createSingleQubitOperationRy(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RY rotation unitary operation**.

The **RY rotation unitary operation** rotates the qubit state around the **Y-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RY rotation unitary operation**.

# See also
- [`ryGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RY rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRy(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912+0.0im  -0.247404+0.0im
0.247404+0.0im   0.968912+0.0im
```
"""
function createSingleQubitOperationRy(theta::Float64)
    # createSingleQubitOperationU3(theta, 0., 0.)
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(pi/2, pi/2)))
end

#

""" 
	rzGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RZ rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RZ rotation gate** performs a unitary operation that rotates the qubit state around the **Z-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RZ rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRz` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRz`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RZ rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rzGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#19#20"{Float64}}(QubiSim.var"#19#20"{Float64}(0.5), [2], "RZ")
```
"""
function rzGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRz(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RZ"
    ))
end

""" 
    createSingleQubitOperationRz(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RZ rotation unitary operation**.

The **RZ rotation unitary operation** rotates the qubit state around the **Z-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RZ rotation unitary operation**.

# See also
- [`rzGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RZ rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRz(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912-0.247404im       0.0+0.0im
0.0+0.0im       0.968912+0.247404im
```
"""
function createSingleQubitOperationRz(theta::Float64)
    # createSingleQubitOperationU3(0., 0., theta) # seems to be wrong
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(0., 0.)))
end

#

""" 
	rotationGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64, 
		phi::Float64, 
		alpha::Float64
	)	-> Nothing

Adds a **rotation gate** around an arbitrary axis to the `quantumCircuit` on the specified `qubit`.

The **rotation gate** around an arbitrary axis performs a unitary operation that rotates the qubit state around the arbitrary axis by an angle ``α``.

The	arbitrary axis lies in the Bloch sphere and is specified by the polar angles ``θ`` and ``φ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRotation` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.
- `phi::Float64` — Rotation angle ``φ`` in radians.
- `alpha::Float64` — Rotation angle ``α`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRotation`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **rotation gate** around the arbitrary axis with polar angles ``θ=0.2``[rad] and ``φ=0.3``[rad] about angle ``α=0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rotationGate!(qc, 2, 0.2, 0.3, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#21#22"{Float64, Float64, Float64}}(QubiSim.var"#21#22"{Float64, Float64, Float64}(0.2, 0.3, 0.5), [2], "Rotation")
```
"""
function rotationGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64, phi::Float64, alpha::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRotation(theta, phi, alpha),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Rotation"
    ))
end

""" 
    createSingleQubitOperationRotation(
		theta::Float64,
		phi::Float64,
		alpha::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **rotation unitary operation** around an arbitrary axis.

The **rotation unitary operation** rotates the qubit state around the arbitrary axis by an angle ``α``.

The	arbitrary axis lies in the Bloch sphere and is specified by the polar angles ``θ`` and ``φ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.
- `phi::Float64` — Rotation angle ``φ`` in radians.
- `alpha::Float64` — Rotation angle ``α`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **rotation unitary operation** around an arbitrary axis.

# See also
- [`rotationGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **rotation unitary operation** around an arbitrary axis with polar angles ``θ=0.2``[rad] and ``φ=0.3``[rad] about angle ``α=0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRotation(0.2, 0.3, 0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
  0.968912-0.242472im   -0.0145253-0.0469563im
 0.0145253-0.0469563im    0.968912+0.242472im
```
"""
function createSingleQubitOperationRotation(theta::Float64, phi::Float64, alpha::Float64)
	UnitaryOperation(exp(-im * (alpha/2.)*createSingleQubitMeasureOperator(theta, phi)))
end

#

""" 
	tGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **T gate** to the `quantumCircuit` on the specified `qubit`.

The **T gate** (also known as the ``π/8`` gate) performs a unitary operation that rotates the phase of the qubit state by ``π/4`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **T unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **T gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> tGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#23#24"}(QubiSim.var"#23#24"(), [2], "T")
```
"""
function tGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationT(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"T"
    ))
end

""" 
    createSingleQubitOperationT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **T unitary operation**.

The **T unitary operation** rotates the phase of the qubit state by ``π/4`` radians about the Z-axis of the Bloch sphere.

This unitary operation is equivalent to applying the ``U_1(π/4)`` or ``U_3(0, 0, π/4)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **T unitary operation**.

# See also
- [`tGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **T unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationT()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.707107+0.707107im
```
"""
function createSingleQubitOperationT()
    createSingleQubitOperationU1(pi/4.)
end

#

""" 
	tdGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **T-dagger gate** to the `quantumCircuit` on the specified `qubit`.

The **T-dagger gate** (also known as the inverse T gate or the ``-π/8`` gate) performs a unitary operation that rotates the phase of the qubit state by ``-π/4`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **T-dagger unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationTd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationTd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **T-dagger gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> tdGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#25#26"}(QubiSim.var"#25#26"(), [2], "Td")
```
"""
function tdGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationTd(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Td"
    ))
end

""" 
    createSingleQubitOperationTd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **T-dagger unitary operation**.

The **T-dagger unitary operation** rotates the phase of the qubit state by ``-π/4`` radians about the Z-axis of the Bloch sphere.

This unitary operation is equivalent to applying the ``U_1(-π/4)`` or ``U_3(0, 0, -π/4)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **T-dagger unitary operation**.

# See also
- [`tdGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **T-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationTd()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.707107-0.707107im
```
"""
function createSingleQubitOperationTd()
    createSingleQubitOperationU1(-pi/4.)
end

#

""" 
	sGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds an **S gate** to the `quantumCircuit` on the specified `qubit`.

The **S gate** (also known as the Phase gate) performs a unitary operation that rotates the phase of the qubit state by ``π/2`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **S unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationS` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationS`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add an **S gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> sGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#27#28"}(QubiSim.var"#27#28"(), [2], "S")
```
"""
function sGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationS(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"S"
    ))
end

""" 
    createSingleQubitOperationS() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **S unitary operation**.

The **S unitary operation** rotates the phase of the qubit state by ``π/2`` radians about the Z-axis of the Bloch sphere.

It is the square root of the **Z unitary operation** (``S² = Z``).

This unitary operation is equivalent to applying the ``U_1(π/2)`` or ``U_3(0, 0, π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **S unitary operation**.

# See also
- [`sGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **S unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationS()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.0+1.0im
```
"""
function createSingleQubitOperationS()
    createSingleQubitOperationU1(pi/2.)
end

#

""" 
	sdGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds an **S-dagger gate** to the `quantumCircuit` on the specified `qubit`.

The **S-dagger gate** (also known as the inverse S or inverse Phase gate) performs a unitary operation that rotates the phase of the qubit state by ``-π/2`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **S-dagger unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationSd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationSd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add an **S-dagger gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> sdGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#29#30"}(QubiSim.var"#29#30"(), [2], "Sd")
```
"""
function sdGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationSd(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Sd"
    ))
end

""" 
    createSingleQubitOperationSd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **S-dagger unitary operation**.

The **S-dagger unitary operation** rotates the phase of the qubit state by ``-π/2`` radians about the Z-axis of the Bloch sphere.

It is the square root of the **Z unitary operation** (``S^†² = Z``).

This unitary operation is equivalent to applying the ``U_1(-π/2)`` or ``U_3(0, 0, -π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **S-dagger unitary operation**.

# See also
- [`sdGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **S-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationSd()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.0-1.0im
```
"""
function createSingleQubitOperationSd()
	createSingleQubitOperationU1(-pi/2.)
end

#

""" 
	projectionGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		v::AbstractVector{<:Number}
	) -> Nothing

Adds a **projection gate** to the `quantumCircuit` on the specified `qubits` that projects onto an arbitrary quantum state vector ``|v⟩``.

The **projection gate** performs a non-unitary operation that projects the system onto the quantum state defined by ``|v⟩``, a complex-valued state vector.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **projection unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationProjection` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` of length ``2^M`` representing the target quantum state for projection.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- `ErrorException` - If the length of ``|v⟩`` does not match ``2^M``.

# Notes
- The vector ``|v⟩`` should be normalized, i.e., ``‖|v⟩‖ ≈ 1`` to ensure the projector behaves correctly.

# See also
- [`createNQubitOperationProjection`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **projection gate** onto the state vector ``|v⟩ = [1; 0; 0; 0]`` to a 2-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> projectionGate!(qc, [1,2], Vector{ComplexF64}([1; 0; 0; 0]))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#31#32"{Vector{ComplexF64}}}(QubiSim.var"#31#32"{Vector{ComplexF64}}(ComplexF64[1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]), [1, 2], "Projection")
```
"""
function projectionGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, v::AbstractVector{<:Number})
	if log2(length(v))==length(qubits)
		addGate!(quantumCircuit, UnitaryGate(
    	    () -> createNQubitOperationProjection(v),
        	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Projection"
    	))
	else
		throw("The length of vector v does not match the specified amount of qubits for the projection gate.")
	end
end

""" 
    createNQubitOperationProjection(
		v::AbstractVector{<:Number}
	) -> UnitaryOperation

Creates a structure encapsulating the matrix representing the **projection operation**.

The **projection operation** projects the system onto the quantum state defined by ``|v⟩``, a complex-valued state vector.

This **projection operation** is not unitary but is returned as a `UnitaryOperation` object for integration into the circuit model.

# Arguments
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` having length equal to a power of two, representing the target quantum state for projection.

# Returns
- `UnitaryOperation` — A structure encapsulating the rank-1 projection matrix ``P = |v⟩⟨v|`` representation of the **projection operation**.

# Notes
- The vector ``|v⟩`` should be normalized, i.e., ``‖|v⟩‖ ≈ 1`` to ensure the projector behaves correctly.
- This operation is not strictly unitary; it is packaged in a UnitaryOperation for compatibility purposes.

# See also
- [`projectionGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the matrix representation of the **projection operation** onto vector ``|v⟩ = [1; 0; 0; 0]``:
```julia-repl
julia> unitaryOperation = createNQubitOperationProjection([1; 0; 0; 0])
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
```
"""
function createNQubitOperationProjection(v::AbstractVector{<:Number})
	UnitaryOperation(v*v')
end

#

""" 
	reflectionGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		v::AbstractVector{<:Number}, 
		theta::Float64
	) -> Nothing

Adds a **reflection gate** to the `quantumCircuit` on the specified `qubit`.

The **reflection gate** performs a unitary operation that reflects the quantum state about the hyperplane orthogonal to an arbitrary quantum state vector ``|v⟩``, scaled by a complex phase factor determined by the angle ``θ``.

This operation implements the generalized reflection operator: ``R(|v⟩, θ) = I - (1 - exp(−i⋅θ))⋅|v⟩⟨v|`` where ``|v⟩`` is a normalized quantum state vector.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **reflection unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationReflection` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` of length ``2^M`` representing the quantum state vector to reflect about.
- `theta::Float64` - The scaling phase angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- `ErrorException` - If the length of ``|v⟩`` is not equal to ``2^M``.

# Notes
- The vector ``|v⟩`` **must be** normalized, i.e., ``‖|v⟩‖ ≈ 1``.

# See also
- [`createNQubitOperationReflection`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **reflection gate** around vector ``|v⟩ = [1; 0; 0; 0]`` with scaling angle ``θ = π`` to a 2-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> reflectionGate!(qc, [1,2], [1;0;0;0], 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#33#34"{Vector{Int64}, Float64}}(QubiSim.var"#33#34"{Vector{Int64}, Float64}([1, 0, 0, 0], 3.141592653589793), [1, 2], "Reflection")
```
"""
function reflectionGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, v::AbstractVector{<:Number}, theta::Float64)
	if log2(length(v))==length(qubits)
		addGate!(quantumCircuit, UnitaryGate(
    	    () -> createNQubitOperationReflection(v, theta),
        	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Reflection"
    	))
	else
		throw("The length of vector v does not match the specified amount of qubits for the reflection gate.")
	end
end

""" 
    createNQubitOperationReflection(
		v::AbstractVector{<:Number}, 
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **reflection unitary operation**.

The **reflection unitary operation** reflects the quantum state about the hyperplane orthogonal to an arbitrary quantum state vector ``|v⟩``, scaled by a complex phase factor determined by the angle ``θ``.

This operation implements the generalized reflection operator: ``R(|v⟩, θ) = I - (1 - exp(−i⋅θ))⋅|v⟩⟨v|`` where ``|v⟩`` is a normalized quantum state vector.

# Arguments
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` having length equal to a power of two, representing the quantum state vector to reflect about.
- `theta::Float64` - The scaling phase angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **reflection unitary operation**.

# Notes
- The vector ``|v⟩`` **must be** normalized, i.e., ``‖|v⟩‖ ≈ 1``.

# See also
- [`reflectionGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **reflection unitary operation** about vector ``|v⟩ = [1; 0; 0; 0]`` with scaling angle ``θ = π``:
```julia-repl
julia> unitaryOperation = createNQubitOperationReflection(ComplexVector([1; 0; 0; 0]), 1.0*pi)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 -1.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im
  0.0-0.0im  1.0-0.0im  0.0-0.0im  0.0-0.0im
  0.0-0.0im  0.0-0.0im  1.0-0.0im  0.0-0.0im
  0.0-0.0im  0.0-0.0im  0.0-0.0im  1.0-0.0im
```
"""
function createNQubitOperationReflection(v::AbstractVector{<:Number}, theta::Float64)
	UnitaryOperation(ComplexMatrix(I, length(v), length(v))-(1-exp(-im*theta))*v*v')
end

#

""" 
	u1Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		lambda::Float64
	) -> Nothing

Adds a **U₁(λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angle ``λ``.

The **U₁(λ) gate** performs a unitary operation that applies a phase rotation by angle ``λ`` around the Z-axis of the Bloch sphere: ``U₁(λ) = U₃(0, 0, λ) = [1, 0; 0, exp(i⋅λ)]``.

This unitary operation is equivalent to `RZ(λ)` up to a global phase.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₁(λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU1` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `lambda::Float64` - The phase rotation angle ``λ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU1`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **``U₁(λ)`` gate** with phase rotation ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u1Gate!(qc, 2, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#35#36"{Float64}}(QubiSim.var"#35#36"{Float64}(3.141592653589793), [2], "U1")
```
"""
function u1Gate!(quantumCircuit::QuantumCircuit, qubit::Int, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createSingleQubitOperationU1(lambda),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U1"
	))
end

""" 
    createSingleQubitOperationU1(
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₁(λ) unitary operation**.

The **U₁(λ) unitary operation** applies a phase rotation by angle ``λ`` around the Z-axis of the Bloch sphere: ``U₁(λ) = U₃(0, 0, λ) = [1, 0; 0, exp(i⋅λ)]``.

This unitary operation is equivalent to `RZ(λ)` up to a global phase.

# Arguments
- `lambda::Float64` - The phase rotation angle ``λ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₁(λ) unitary operation**.

# See also
- [`u1Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₁(λ) unitary operation** with phase rotation ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU1(1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0-0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function createSingleQubitOperationU1(lambda::Float64)
    createSingleQubitOperationU(0., 0., lambda)
end

#

""" 
	u2Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		phi::Float64, 
		lambda::Float64
	) -> Nothing

Adds a **U₂(φ, λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angles ``φ`` and ``λ``.

The **U₂(φ, λ) gate** performs a unitary operation that applies a double phase rotation defined by: ``U₂(φ, λ) = U₃(π/2, φ, λ) = 1/√2⋅[1, -exp(i⋅λ); exp(i⋅φ), exp(i⋅(φ+λ))]``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₂(φ, λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU2` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `phi::Float64` -  The first phase rotation angle ``φ`` in radians.
- `lambda::Float64` - The second phase rotation angle ``λ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU2`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **U₂(φ, λ) gate** with phase rotations ``φ = 0.5π`` and ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u2Gate!(qc, 2, 0.5*pi, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#37#38"{Float64, Float64}}(QubiSim.var"#37#38"{Float64, Float64}(1.5707963267948966, 3.141592653589793), [2], "U2")
```
"""
function u2Gate!(quantumCircuit::QuantumCircuit, qubit::Int, phi::Float64, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createSingleQubitOperationU2(phi, lambda),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U2"
	))
end

""" 
    createSingleQubitOperationU2(
		phi::Float64, 
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₂(φ, λ) unitary operation**.

The **U₂(φ, λ) unitary operation** applies a double phase rotation defined by: ``U₂(φ, λ) = U₃(π/2, φ, λ) = 1/√2⋅[1, -exp(i⋅λ); exp(i⋅φ), exp(i⋅(φ+λ))]``.

# Arguments
- `phi::Float64` - The first phase rotation angle ``φ`` in radians.
- `lambda::Float64` - The second phase rotation angle ``λ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₂(φ, λ) unitary operation**.

# See also
- [`u2Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₂(φ, λ) unitary operation** with phase rotations ``φ = 0.5π`` and ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU2(0.5*pi, 1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
    0.707107+0.0im  0.707107+0.0im
    0.0+0.707107im  0.0-0.707107im
```
"""
function createSingleQubitOperationU2(phi::Float64, lambda::Float64)
    createSingleQubitOperationU(pi/2,phi,lambda)
end

#

""" 
	u3Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64, 
		phi::Float64, 
		lambda::Float64
	) -> Nothing

Adds a **U₃(θ, φ, λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angles ``θ``, ``φ`` and ``λ``.

The **U₃(θ, φ, λ) gate** performs the most general unitary operation that applies a triple phase rotation defined by: ``U₃(θ, φ, λ) = [cos(θ/2), -exp(i⋅λ)⋅sin(θ/2); exp(i⋅φ)⋅sin(θ/2), exp(i⋅(φ+λ))⋅cos(θ/2)]``

Special cases:
- ``U₃(0, 0, 0)`` = Identity
- ``U₃(π/2, φ, λ) = U₂(φ, λ)``
- ``U₃(0, 0, λ) = U₁(λ)``

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₃(θ, φ, λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU3` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` - Rotation angle ``θ`` around the Y-axis in radians.
- `phi::Float64` - Phase angle ``φ`` for pre-Z rotation in radians.
- `lambda::Float64` - Phase angle ``λ`` for post-Z rotation in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU3`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **U₃(θ, φ, λ) gate** with phase rotations ``θ = 0.25π``, ``φ = 0.5π`` and ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u3Gate!(qc, 2, 0.25*pi, 0.5*pi, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#39#40"{Float64, Float64, Float64}}(QubiSim.var"#39#40"{Float64, Float64, Float64}(0.7853981633974483, 1.5707963267948966, 3.141592653589793), [2], "U3")
```
"""
function u3Gate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64, phi::Float64, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationU3(theta, phi, lambda),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U3"
    ))
end

""" 
    createSingleQubitOperationU3(
		theta::Float64, 
		phi::Float64, 
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₃(θ, φ, λ) unitary operation**.

The **U₃(θ, φ, λ) unitary operation** is the most general unitary operation that applies a triple phase rotation defined by: ``U₃(θ, φ, λ) = [cos(θ/2), -exp(i⋅λ)⋅sin(θ/2); exp(i⋅φ)⋅sin(θ/2), exp(i⋅(φ+λ))⋅cos(θ/2)]``

Special cases:
- ``U₃(0, 0, 0)`` = Identity
- ``U₃(π/2, φ, λ) = U₂(φ, λ)``
- ``U₃(0, 0, λ) = U₁(λ)``

# Arguments
- `theta::Float64` - Rotation angle ``θ`` around the Y-axis in radians.
- `phi::Float64` - Phase angle ``φ`` for pre-Z rotation in radians.
- `lambda::Float64` - Phase angle ``λ`` for post-Z rotation in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₃(θ, φ, λ) unitary operation**.

# See also
- [`u3Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₃(θ, φ, λ) unitary operation** with phase rotations ``θ = 0.25π``, ``φ = 0.5π`` and ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU3(0.25*pi, 0.5*pi, 1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
    0.92388+0.0im   0.382683+0.0im
    0.0+0.382683im  0.0-0.92388im
```
"""
function createSingleQubitOperationU3(theta::Float64, phi::Float64, lambda::Float64)
    createSingleQubitOperationU(theta, phi, lambda)
end

#

""" 
	cnotGate!(
		quantumCircuit::QuantumCircuit, 
		cqubit::Int, 
		tqubit::Int
	) -> Nothing

Adds a **CNOT gate** to the `quantumCircuit` on the specified control qubit `cqubit` and target qubit `tqubit`.

The **CNOT gate** (also known as the Controlled-NOT gate) performs a unitary operation that flips the target qubit (``|0⟩ ↔ |1⟩``) *if and only if* the control qubit is in the state ``|1⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **CNOT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationCNOT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `cqubit::Int` — Index of the control qubit (interpreted using the circuit’s `settings`).
- `tqubit::Int` — Index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationCNOT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **CNOT gate** to a 3-qubit quantum circuit on control qubit 1 and target qubit 2:
# Example
Add a CNOT gate  to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> cnotGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#41#42"}(QubiSim.var"#41#42"(), [1, 2], "CNOT")
```
"""
function cnotGate!(quantumCircuit::QuantumCircuit, cqubit::Int, tqubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationCNOT(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [cqubit,tqubit]),
		"CNOT"
	))
end

""" 
    createDoubleQubitOperationCNOT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **CNOT unitary operation**.

The **CNOT unitary operation** flips the target qubit (``|0⟩ ↔ |1⟩``) *if and only if* the control qubit is in the state ``|1⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **CNOT unitary operation**.

# See also
- [`cnotGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **CNOT unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationCNOT()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationCNOT()
	UnitaryOperation([1. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 1.;0. 0. 1. 0.])
end

#

""" 
	cnotReverseGate!(
		quantumCircuit::QuantumCircuit, 
		cqubit::Int, 
		tqubit::Int
	) -> Nothing

Adds a **reversed CNOT gate** to the `quantumCircuit` on the specified control qubit `cqubit` and target qubit `tqubit`.

The **reversed CNOT gate** (also known as a target-controlled NOT gate) performs a unitary operation that flips the control qubit (``|0⟩ ↔ |1⟩``) *if and only if* the target qubit is in the state ``|1⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **reversed CNOT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationCNOTReverse` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `cqubit::Int` — Index of the control qubit (interpreted using the circuit’s `settings`).
- `tqubit::Int` — Index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationCNOTReverse`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **reversed CNOT gate** to a 3-qubit quantum circuit on control qubit 1 and target qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> cnotReverseGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#43#44"}(QubiSim.var"#43#44"(), [1, 2], "CNOT-Reverse")
```
"""
function cnotReverseGate!(quantumCircuit::QuantumCircuit, cqubit::Int, tqubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationCNOTReverse(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [cqubit,tqubit]),
		"CNOT-Reverse"
	))
end

""" 
    createDoubleQubitOperationCNOTReverse() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **reversed CNOT unitary operation**.

The **reversed CNOT unitary operation** flips the control qubit (``|0⟩ ↔ |1⟩``) *if and only if* the target qubit is in the state ``|1⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **reversed CNOT unitary operation**.

# See also
- [`cnotReverseGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **reversed CNOT unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationCNOTReverse()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationCNOTReverse()
    UnitaryOperation([1. 0. 0. 0.;0. 0. 0. 1.;0. 0. 1. 0.;0. 1. 0. 0.])
end

#

""" 
	swapGate!(
		quantumCircuit::QuantumCircuit, 
		qubit1::Int, 
		qubit2::Int
	) -> Nothing

Adds a **SWAP gate** to the `quantumCircuit` on the specified qubits `qubit1` and `qubit2`.

The **SWAP gate** performs a unitary operation that fully exchanges the basis states between two qubits: ``|a⟩⊗|b⟩ → |b⟩⊗|a⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **SWAP unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationSWAP` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit1::Int` — Index of the first qubit (interpreted using the circuit’s `settings`).
- `qubit2::Int` — Index of the second qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationSWAP`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **SWAP gate** to a 3-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> swapGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")
```
"""
function swapGate!(quantumCircuit::QuantumCircuit, qubit1::Int, qubit2::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationSWAP(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit1,qubit2]),
		"SWAP"
	))
end

""" 
    createDoubleQubitOperationSWAP() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **SWAP unitary operation**.

The **SWAP unitary operation** fully exchanges the basis states between two qubits: ``|a⟩⊗|b⟩ → |b⟩⊗|a⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **SWAP unitary operation**.

# See also
- [`swapGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **SWAP unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationSWAP()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
```
"""
function createDoubleQubitOperationSWAP()
    UnitaryOperation([1. 0. 0. 0.;0. 0. 1. 0.;0. 1. 0. 0.;0. 0. 0. 1.])
end

#

""" 
	phaseGate!(
		quantumCircuit::QuantumCircuit, 
		qubit1::Int, 
		qubit2::Int, 
		phi::Float64
	) -> Nothing

Adds a **controlled phase gate** to the `quantumCircuit` on the specified qubits `qubit1` and `qubit2`, with a phase angle ``ϕ``.

The **controlled phase gate** performs a unitary operation that applies a conditional phase rotation ``|11⟩ ↦ exp(i⋅ϕ)⋅|11⟩`` while leaving the other computational basis states unchanged.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled phase unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationPhase` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit1::Int` — Index of the first qubit (interpreted using the circuit’s `settings`).
- `qubit2::Int` — Index of the second qubit (interpreted using the circuit’s `settings`).
- `phi::Float64` — Phase angle ``ϕ`` in radians (i.e., the amount of phase applied to the ``|11⟩`` state).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationPhase`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **controlled phase gate** with phase angle ``ϕ = 0.5π`` on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> phaseGate!(qc, 1, 2, 0.5*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(1.5707963267948966), [1, 2], "Phase")
```
"""
function phaseGate!(quantumCircuit::QuantumCircuit, qubit1::Int, qubit2::Int, phi::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationPhase(phi),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit1,qubit2]),
		"Phase"
	))
end

""" 
    createDoubleQubitOperationPhase(
		phi::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled phase unitary operation**.

The **controlled phase unitary operation** applies a conditional phase rotation ``|11⟩ ↦ exp(i⋅ϕ)⋅|11⟩`` while leaving the other computational basis states unchanged.

# Arguments
- `phi::Float64` — Phase angle ``ϕ`` in radians (i.e., the amount of phase applied to the ``|11⟩`` state).

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled phase unitary operation**.

# See also
- [`phaseGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **controlled phase unitary operation** with phase angle ``ϕ = 0.5π``:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationPhase(0.5*pi)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+1.0im
```
"""
function createDoubleQubitOperationPhase(phi::Float64)
	UnitaryOperation([1. 0. 0. 0.;0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. exp(im*phi)])
end

#

""" 
	unitaryUGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		U::UnitaryOperation; 
		accuracyCheckForUnitarity = 10*eps(1.0)
	) -> Nothing

Adds a **unitary gate** to the `quantumCircuit` on the specified `qubits` with a user-supplied **unitary operation** ``U``.

The **unitary gate** performs the user-supplied **unitary operation** ``U``.

The function performs an internal unitarity check to ensure that the user-supplied **unitary operation** ``U`` satisfies the condition ``U^†⋅U ≈ I`` within a specified numerical tolerance.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **unitary operation**, and adds 
it to the circuit. During compilation, this factory function simply returns the user-suppied **unitary operation** `U` as the produced unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices (interpreted using the circuit’s `settings`) on which the gate is applied.
- `U::UnitaryOperation` — A structure encapsulating the unitary matrix ``U`` representation of the **unitary operation**.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for unitarity validation, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if ``U`` is not unitary within the given tolerance.

# See also
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Supported **unitary operations** include
- Single-qubit gates:
  - [`createSingleQubitOperationU1`](@ref)
  - [`createSingleQubitOperationU2`](@ref)
  - [`createSingleQubitOperationU3`](@ref)
  - [`createSingleQubitOperationX`](@ref)
  - [`createSingleQubitOperationY`](@ref)
  - [`createSingleQubitOperationZ`](@ref)
  - [`createSingleQubitOperationH`](@ref)
  - [`createSingleQubitOperationId`](@ref)
  - [`createSingleQubitOperationRx`](@ref)
  - [`createSingleQubitOperationRy`](@ref)
  - [`createSingleQubitOperationRz`](@ref)
  - [`createSingleQubitOperationRotation`](@ref)
  - [`createSingleQubitOperationT`](@ref)
  - [`createSingleQubitOperationTd`](@ref)
  - [`createSingleQubitOperationS`](@ref)
  - [`createSingleQubitOperationSd`](@ref)
- Double-qubit gates:
  - [`createDoubleQubitOperationCNOT`](@ref)
  - [`createDoubleQubitOperationCNOTReverse`](@ref)
  - [`createDoubleQubitOperationSWAP`](@ref)
  - [`createDoubleQubitOperationPhase`](@ref)
- Multi-qubit gates:
  - [`createNQubitOperationProjection`](@ref)
  - [`createNQubitOperationReflection`](@ref)
  - [`createNQubitOperationExpH`](@ref)
  - [`createNQubitOperationQFT`](@ref)
  - [`createNQubitOperationIQFT`](@ref)
  - [`createNQubitOperationControlledU`](@ref)
  - [`createNQubitOperationQPE`](@ref)
  - [`compileToSingleGate`](@ref)

# Example
Add a custom **unitary gate** with a **unitary operation** for a Z-gate on qubit 1 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> unitaryUGate!(qc, [1], createSingleQubitOperationZ())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#50#51"{UnitaryOperation}}(QubiSim.var"#50#51"{UnitaryOperation}(UnitaryOperation(ComplexF64[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 1.2246467991473532e-16im])), [1], "Unitary")
```
"""
function unitaryUGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, U::UnitaryOperation; accuracyCheckForUnitarity = 10*eps(1.0))
	if all(abs.(U.U'*U.U-ComplexMatrix(I,size(U.U))).<accuracyCheckForUnitarity)
		addGate!(quantumCircuit, UnitaryGate(
		    () -> U,
    		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Unitary"
		))
	else
		throw("U:$U is not unitary")
	end
end

#

""" 
	expHGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		H::ComplexMatrix, 
		t_over_h::Float64; 
		accuracyCheckForHermiticity = 10*eps(1.0))
	) -> Nothing

Adds a **time evolution gate** to the `quantumCircuit` on the specified `qubits` for a Hamiltonian ``H`` and scaled time ``t/ℏ``.

The **time evolution gate** performs a unitary operation that applies the time evolution:
```math
exp(-i·H·t/ℏ)
```

where ``H`` is the system's Hamiltonian, and ``t/ℏ`` is the ratio of evolution time ``t`` over Planck’s constant ``ℏ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **time evolution unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationExpH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `H::ComplexMatrix` — The Hamiltonian matrix ``H`` of size ``2^M × 2^M`` governing the evolution, must be Hermitian.
- `t_over_h::Float64` — Time evolution parameter ``t/ℏ``, expressed in radians.
- `accuracyCheckForHermiticity::Float64` (optional) — Numerical tolerance for Hermiticity validation of ``H``, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of ``H`` are incompatible with the number of target qubits or fails the Hermiticity check (i.e., ``H^† ≠ H`` within the specified tolerance).

# See also
- [`createNQubitOperationExpH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **time evolution gate** for the Hamiltonian:

    H = [ 0.460075592255305    -1.110720734539591;
         -1.110720734539591     2.681517061334488 ]

with evolution time ``t/ℏ = 1.0`` on qubit 1 to a 2-qubit circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> expHGate!(qc, [1], ComplexMatrix([0.460075592255305 -1.110720734539591;-1.110720734539591 2.681517061334488]), 1.0)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#53#54"{Matrix{ComplexF64}, Float64}}(QubiSim.var"#53#54"{Matrix{ComplexF64}, Float64}(ComplexF64[0.460075592255305 + 0.0im -1.110720734539591 + 0.0im; -1.110720734539591 + 0.0im 2.681517061334488 + 0.0im], 1.0), [1], "expH")
```
"""
function expHGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, H::ComplexMatrix, t_over_h::Float64; accuracyCheckForHermiticity = 10*eps(1.0))
	if size(H) == (2^length(qubits), 2^length(qubits))
		if all(abs.(H'-H).<accuracyCheckForHermiticity)
			addGate!(quantumCircuit, UnitaryGate(
	    		() -> createNQubitOperationExpH(H, t_over_h),
    			convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
				"expH"
			))
		else
			throw("H is not Hermitian")
		end
	else
		throw("Dimensions of H are incompatible with the number of target qubits")
	end
end

""" 
	createNQubitOperationExpH(
		H::ComplexMatrix, 
		t_over_h::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **time evolution unitary operation**.

The **time evolution unitary operation** applies the time evolution: 
```math
exp(-i·H·t/ℏ)
```

where ``H`` is the system's Hamiltonian, and ``t/ℏ`` is the ratio of evolution time ``t`` over Planck’s constant ``ℏ``.

# Arguments
- `H::ComplexMatrix` — The Hamiltonian matrix ``H`` governing the evolution, must be Hermitian.
- `t_over_h::Float64` — Time evolution parameter ``t/ℏ``, expressed in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **time evolution unitary operation**.

# Notes
- This function supports any ``M``-qubit system as long as the size of ``H`` is ``2^M × 2^M``.
- ``H`` must be Hermitian to ensure the resulting operation is unitary.

# See also
- [`expHGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **time evolution unitary operation** for a single-qubit system with the Hamiltonian:

    H = [ 0.460075592255305    -1.110720734539591;
         -1.110720734539591     2.681517061334488 ]

and evolution time ``t/ℏ = 1.0``:
```julia-repl
julia> unitaryOperation = createNQubitOperationExpH(ComplexMatrix([0.460075592255305  -1.110720734539591;-1.110720734539591   2.681517061334488]), 1.0)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
function createNQubitOperationExpH(H::ComplexMatrix,t_over_h::Float64)
	UnitaryOperation(exp(-im * H * t_over_h))
end

#

""" 
	qftGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits
	) -> Nothing

Adds a **Quantum Fourier Transform (QFT) gate** to the `quantumCircuit` on the specified `qubits`.

The **QFT gate** performs a unitary operation that transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to the length of `qubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **QFT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationQFT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createNQubitOperationQFT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **QFT gate** on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> qftGate!(qc, [1,2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#55#56"{Vector{Int64}}}(QubiSim.var"#55#56"{Vector{Int64}}([1, 2]), [1, 2], "QFT")
```
"""
function qftGate!(quantumCircuit::QuantumCircuit, qubits::Qubits)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createNQubitOperationQFT(length(qubits)),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
		"QFT"
	))
end

""" 
	createNQubitOperationQFT(
		numberOfQubits::Int
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Quantum Fourier Transform (QFT) unitary operation** on `numberOfQubits` qubits.

The **QFT unitary operation** transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

Under the hood, it builds the full quantum circuit implementing the **QFT** and 
compiles it into a single unitary matrix representing the **QFT unitary operation**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **QFT** is applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **QFT unitary operation**.

# See also
- [`qftGate!`](@ref)
- [`createNQubitQFTQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **QFT unitary operation** acting on 2 qubits:
```julia-repl
julia> unitaryOperation = createNQubitOperationQFT(2)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 0.5+0.0im   0.5+0.0im   0.5+0.0im   0.5+0.0im
 0.5+0.0im   0.0+0.5im  -0.5+0.0im   0.0-0.5im
 0.5+0.0im  -0.5+0.0im   0.5+0.0im  -0.5+0.0im
 0.5+0.0im   0.0-0.5im  -0.5+0.0im   0.0+0.5im
```
"""
function createNQubitOperationQFT(numberOfQubits::Int)
	quantumCircuit = createNQubitQFTQuantumCircuit(numberOfQubits)
	return compileToSingleGate(quantumCircuit)
end

"""
    createNQubitQFTQuantumCircuit(
		numberOfQubits::Int
	) -> QuantumCircuit{Z, I}

Builds the quantum circuit implementing the **Quantum Fourier Transform (QFT)** on `numberOfQubits` qubits.

The **QFT** transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **QFT** is applied.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **QFT**

# See also
- [`createNQubitOperationQFT`](@ref)
- [`QuantumCircuit`](@ref)

# Example
Create the quantum circuit of a 2-qubit **QFT**:
```julia-repl
julia> quantumCircuit = createNQubitQFTQuantumCircuit(2)
julia> quantumCircuit.circuit
4-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(1.5707963267948966), [2, 1], "Phase")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")])
```
"""
function createNQubitQFTQuantumCircuit(numberOfQubits::Int)
	quantumCircuit = createQuantumCircuit(numberOfQubits)
	
	for m in 1:numberOfQubits
        hGate!(quantumCircuit, m)
		for n in 1:(numberOfQubits - m)
            phaseGate!(quantumCircuit, n + m, m, pi / (2 ^ n))
		end
	end
	for m in 1:Int(floor(numberOfQubits / 2))
        swapGate!(quantumCircuit, m, numberOfQubits - m + 1)
	end
	return quantumCircuit
end

#

""" 
	iqftGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits
	) -> Nothing

Adds an **Inverse Quantum Fourier Transform (IQFT) gate** to the `quantumCircuit` on the specified `qubits`.

The **IQFT gate** performs a unitary operation that transforms the basis state ``|x⟩`` into: 
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to the length of `qubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **IQFT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationIQFT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createNQubitOperationIQFT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add an **IQFT gate** on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc=createQuantumCircuit(2)
julia> iqftGate!(qc, [1,2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#57#58"{Vector{Int64}}}(QubiSim.var"#57#58"{Vector{Int64}}([1, 2]), [1, 2], "IQFT")
```
"""
function iqftGate!(quantumCircuit::QuantumCircuit, qubits::Qubits)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createNQubitOperationIQFT(length(qubits)),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
		"IQFT"
	))
end

""" 
    createNQubitOperationIQFT(
		numberOfQubits::Int
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Inverse Quantum Fourier Transform (IQFT) unitary operation** on `numberOfQubits` qubits.

The **IQFT unitary operation** transforms the basis state ``|x⟩`` into:
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

Under the hood, it builds the full quantum circuit implementing the **IQFT** and 
compiles it into a single unitary matrix representing the **IQFT unitary operation**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **IQFT** is applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **IQFT unitary operation**.

# See also
- [`iqftGate!`](@ref)
- [`createNQubitIQFTQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **IQFT unitary operation** acting on 2 qubits:
```julia-repl
julia> unitaryOperation=createNQubitOperationIQFT(2)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 0.5+0.0im   0.5+0.0im   0.5+0.0im   0.5+0.0im
 0.5+0.0im   0.0-0.5im  -0.5+0.0im   0.0+0.5im
 0.5+0.0im  -0.5+0.0im   0.5+0.0im  -0.5+0.0im
 0.5+0.0im   0.0+0.5im  -0.5+0.0im   0.0-0.5im
```
"""
function createNQubitOperationIQFT(numberOfQubits::Int)
	quantumCircuit = createNQubitIQFTQuantumCircuit(numberOfQubits)
	return compileToSingleGate(quantumCircuit)
end

""" 
    createNQubitIQFTQuantumCircuit(
		numberOfQubits::Int
	) -> QuantumCircuit{Z, I}

Builds the quantum circuit implementing the **Inverse Quantum Fourier Transform (IQFT)** on `numberOfQubits` qubits.

The **IQFT** transforms the basis state ``|x⟩`` into: 
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **IQFT** is applied.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **IQFT**

# See also
- [`createNQubitOperationIQFT`](@ref)
- [`QuantumCircuit`](@ref)

# Example
Create the quantum circuit of a 2-qubit **IQFT**:
```julia-repl
julia> quantumCircuit = createNQubitIQFTQuantumCircuit(2)
julia> quantumCircuit.circuit
4-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(-1.5707963267948966), [2, 1], "Phase")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
```
"""
function createNQubitIQFTQuantumCircuit(numberOfQubits::Int)
	quantumCircuit = createQuantumCircuit(numberOfQubits)

	for m in 1:Int(floor(numberOfQubits / 2))
       swapGate!(quantumCircuit, m, numberOfQubits - m + 1)
	end
	for m in numberOfQubits:-1:1
		for n in 1:(numberOfQubits - m)
            phaseGate!(quantumCircuit, n + m, m, -pi / (2 ^ n))
		end
        hGate!(quantumCircuit, m)
	end
	return quantumCircuit
end

#

""" 
	controlledUGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubits::Qubits, 
		targetQubits::Qubits, 
		UTarget::UnitaryOperation
	) -> Nothing

Adds a **controlled unitary gate** to the `quantumCircuit` on the specified qubits `controlQubits` and `targetQubits` and unitary operation `UTarget`.

The **controlled unitary gate** performs a unitary operation that applies the unitary operation `UTarget` on the `targetQubits` if the `controlQubits` are **all** in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationControlledU` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubits::Qubits` — Control qubit indices (interpreted via the circuit's `settings`).
- `targetQubits::Qubits` — Target qubit indices (interpreted via the circuit's `settings`).
- `UTarget::UnitaryOperation` — The unitary operation to be applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of `UTarget` are incompatible with the number of `targetQubits`.

# See also
- [`createNQubitOperationControlledU`](@ref)
- [`UnitaryOperation`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **controlled Pauli-Z gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledUGate!(qc, [1], [2], createSingleQubitOperationZ())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [2], UnitaryOperation(ComplexF64[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 1.2246467991473532e-16im])), [1, 2], "Controlled-Unitary")
```
"""
function controlledUGate!(quantumCircuit::QuantumCircuit, controlQubits::Qubits, targetQubits::Qubits, UTarget::UnitaryOperation)
	if size(UTarget.U) == (2^length(targetQubits), 2^length(targetQubits))
		addGate!(quantumCircuit, UnitaryGate(
			() -> createNQubitOperationControlledU(length(controlQubits), length(targetQubits), UTarget),
			convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, union(controlQubits, targetQubits)),
			"Controlled-Unitary"
		))
	else
		throw("Dimensions of UTarget are incompatible with the number of target qubits")
	end
end

""" 
    createNQubitOperationControlledU(
		numberOfControlQubits::Int, 
		numberOfTargetQubits::Int, 
		UTarget::UnitaryOperation
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled unitary operation**.

The **controlled unitary operation** applies the unitary operation `UTarget` on the target qubits if the control qubits are **all** in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubits.

# Arguments
- `numberOfControlQubits::Int` — Number of control qubits.
- `numberOfTargetQubits::Int` — Number of target qubits.
- `UTarget::UnitaryOperation` — The unitary operation to be applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled unitary operation**.

# See also
- [`controlledUGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Z unitary operation** with one control and one target qubit:
```julia-repl
julia> unitaryOperation = createNQubitOperationControlledU(1, 1, createSingleQubitOperationZ())
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im   0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  -1.0+0.0im
```
"""
function createNQubitOperationControlledU(numberOfControlQubits::Int, numberOfTargetQubits::Int, UTarget::UnitaryOperation)
	totalNumberOfQubits = numberOfControlQubits+numberOfTargetQubits
	controlledU = createNQubitOperationId(totalNumberOfQubits)
    controlledU.U[2^totalNumberOfQubits-size(UTarget.U,1)+1:end, 
                  2^totalNumberOfQubits-size(UTarget.U,2)+1:end] = UTarget.U

	return controlledU
end

#

""" 
	controlledHGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Hadamard gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Hadamard gate** performs a unitary operation that applies the Hadamard transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Hadamard operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`hGate!`](@ref)

# Example
Add a **controlled Hadamard gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledHGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#61#62"}(QubiSim.var"#61#62"(), [1, 2], "Controlled-H")
```
"""
function controlledHGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledH(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-H"
	))
end

""" 
    createDoubleQubitOperationControlledH() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Hadamard operation**.

The **controlled Hadamard operation** applies the Hadamard transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Hadamard operation**.

# See also
- [`controlledHGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationH`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Hadamard unitary operation**:
julia> unitaryOperation = createDoubleQubitOperationControlledH()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im       0.0+0.0im        0.0+0.0im
 0.0+0.0im  1.0+0.0im       0.0+0.0im        0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.707107+0.0im   0.707107+0.0im
 0.0+0.0im  0.0+0.0im  0.707107+0.0im  -0.707107+0.0im
```
"""
function createDoubleQubitOperationControlledH()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationH()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledXGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-X gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-X gate** performs a unitary operation that applies the Pauli-X transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-X operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledX` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledX`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`xGate!`](@ref)

# Example
Add a **controlled Pauli-X gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledXGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#63#64"}(QubiSim.var"#63#64"(), [1, 2], "Controlled-X")
```
"""
function controlledXGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledX(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-X"
	))
end

""" 
    createDoubleQubitOperationControlledX() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-X operation**.

The **controlled Pauli-X operation** applies the Pauli-X transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-X operation**.

# See also
- [`controlledXGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationX`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-X unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledX()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationControlledX()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationX()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledYGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-Y gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-Y gate** performs a unitary operation that applies the Pauli-Y transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-Y operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledY` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledY`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`hGate!`](@ref)

# Example
Add a **controlled Pauli-Y gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledYGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#65#66"}(QubiSim.var"#65#66"(), [1, 2], "Controlled-Y")
```
"""
function controlledYGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledY(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Y"
	))
end

""" 
    createDoubleQubitOperationControlledY() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-Y operation**.

The **controlled Pauli-Y operation** applies the Pauli-Y transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-Y operation**.

# See also
- [`controlledYGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationY`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Y unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledY()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0-1.0im
 0.0+0.0im  0.0+0.0im  0.0+1.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationControlledY()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationY()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledZGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-Z gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-Z gate** performs a unitary operation that applies the Pauli-Z transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-Z operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledZ` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledZ`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`zGate!`](@ref)

# Example
Add a **controlled Pauli-Z gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledZGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#67#68"}(QubiSim.var"#67#68"(), [1, 2], "Controlled-Z")
```
"""
function controlledZGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledZ(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Z"
	))
end

""" 
    createDoubleQubitOperationControlledZ() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-Z operation**.

The **controlled Pauli-Z operation** applies the Pauli-Z transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-Z operation**.

# See also
- [`controlledZGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationZ`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Z unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledZ()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im   0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  -1.0+0.0im
```
"""
function createDoubleQubitOperationControlledZ()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationZ()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledTGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled T gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled T gate** performs a unitary operation that applies the T transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled T operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`tGate!`](@ref)

# Example
Add a **controlled T gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledTGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#69#70"}(QubiSim.var"#69#70"(), [1, 2], "Controlled-T")
```
"""
function controlledTGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledT(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-T"
	))
end

""" 
    createDoubleQubitOperationControlledT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled T operation**.

The **controlled T operation** applies the T transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled T operation**.

# See also
- [`controlledTGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationT`](@ref)

# Example
Extract the unitary matrix representation of the **controlled T unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledT()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im      -0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.707107+0.707107im
```
"""
function createDoubleQubitOperationControlledT()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationT()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledTdGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled T-dagger gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled T-dagger gate** performs a unitary operation that applies the T-dagger transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled T-dagger operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledTd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledTd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`tdGate!`](@ref)

# Example
Add a **controlled T-dagger gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledTdGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#71#72"}(QubiSim.var"#71#72"(), [1, 2], "Controlled-Td")
```
"""
function controlledTdGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledTd(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Td"
	))
end

""" 
    createDoubleQubitOperationControlledTd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled T-dagger operation**.

The **controlled T-dagger operation** applies the T-dagger transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled T-dagger operation**.

# See also
- [`controlledTdGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationTd`](@ref)

# Example
Extract the unitary matrix representation of the **controlled T-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledTd()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im      -0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.707107-0.707107im
```
"""
function createDoubleQubitOperationControlledTd()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationTd()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledSGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled S gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled S gate** performs a unitary operation that applies the S transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled S operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledS` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledS`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`sGate!`](@ref)

# Example
Add a **controlled S gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledSGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#73#74"}(QubiSim.var"#73#74"(), [1, 2], "Controlled-S")
```
"""
function controlledSGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledS(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-S"
	))
end

""" 
    createDoubleQubitOperationControlledS() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled S operation**.

The **controlled S operation** applies the S transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled S operation**.

# See also
- [`controlledSGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationS`](@ref)

# Example
Extract the unitary matrix representation of the **controlled S unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledS()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  -0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+1.0im
```
"""
function createDoubleQubitOperationControlledS()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationS()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledSdGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled S-dagger gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled S-dagger gate** performs a unitary operation that applies the S-dagger transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled S-dagger operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledSd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledSd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`sdGate!`](@ref)

# Example
Add a **controlled S-dagger gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledSdGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#75#76"}(QubiSim.var"#75#76"(), [1, 2], "Controlled-Sd")
```
"""
function controlledSdGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledSd(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Sd"
	))
end

""" 
    createDoubleQubitOperationControlledSd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled S-dagger operation**.

The **controlled S-dagger operation** applies the S-dagger transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled S-dagger operation**.

# See also
- [`controlledSdGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationSd`](@ref)

# Example
Extract the unitary matrix representation of the **controlled S-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledSd()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0-1.0im
```
"""
function createDoubleQubitOperationControlledSd()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationSd()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	qpeGate!(
		quantumCircuit::QuantumCircuit, 
		estimationQubits::Qubits,
		targetQubits::Qubits,
		U::UnitaryOperation;
		accuracyCheckForUnitarity = 10 * eps(1.0)
	) -> Nothing

Adds a **Quantum Phase Estimation (QPE) gate** to `quantumCircuit`, acting on the specified `estimationQubits` and `targetQubits`, for the given unitary operation ``U``.

The **QPE gate** estimates the eigenphase(s) of the unitary operation ``U``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **QPE unitary operation**, and adds
it to the circuit. During compilation, this factory function calls `createNQubitOperationQPE` to produce the unitary matrix.

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (`estimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (`targetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on a number of qubits equal to the length of `targetQubits`.
- A numerical check ensures that ``U`` satisfies unitarity: ``U^†·U ≈ I`` within the specified numerical tolerance.
- If the matrix dimension of ``U`` does not match the number of `targetQubits`, or ``U`` fails the unitarity check, an exception is thrown.
- The estimation and target qubit indexing is interpreted using the circuit’s `settings`.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `estimationQubits::Qubits` — Indices of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `targetQubits::Qubits` — Indices of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for verifying the unitarity of ``U``, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of ``U`` are incompatible with the number of `targetQubits` or fails the unitarity check.

# See also
- [`createNQubitOperationQPE`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`UnitaryOperation`](@ref)
- [`Qubits`](@ref)
- [`Settings`](@ref)

# Example
Use a **QPE gate** to estimate the eigenphase of a **phase gate** with 4-qubit accuracy:
```julia-repl
julia> accuracy = 4
julia> phasePart = 11/15
julia> theta = 2*pi*phasePart
julia> qc = createQuantumCircuit(accuracy + 1)
julia> xGate!(qc, accuracy + 1) # to create eigenstate |1> with eigenphase theta from initial state |0>
julia> qpeGate!(qc, Vector(1:accuracy), [accuracy + 1], createSingleQubitOperationU1(theta))
julia> measureGate!(qc, Vector(1:accuracy))
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
julia> iqs = createInitialQubitState(vector, repeat([1. 0.], outer=[accuracy + 1, 1]))
julia> qo = runQuantumProgram(qp, iqs, 1000)
julia> phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^accuracy)
julia> probeMeasureOutcome(qo, 3, "Estimated phase of a "*string(round(theta, digits=3))*" [rad] phase gate: "*string(round(phaseEstimate, digits=3))*" [rad]")
```
"""
function qpeGate!(quantumCircuit::QuantumCircuit, estimationQubits::Qubits, targetQubits::Qubits, U::UnitaryOperation; accuracyCheckForUnitarity = 10*eps(1.0))
	if size(U.U) == (2^length(targetQubits), 2^length(targetQubits))
		if all(abs.(U.U'*U.U-ComplexMatrix(I,size(U.U))).<accuracyCheckForUnitarity)
			addGate!(quantumCircuit, UnitaryGate(
				() -> createNQubitOperationQPE(length(estimationQubits), length(targetQubits), U, quantumCircuit.settings),
				convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, union(estimationQubits, targetQubits)),
				"QPE"
			))
		else
			throw("U:$U is not unitary")
		end
	else
		throw("Dimensions of U are incompatible with the number of target qubits")
	end
end

""" 
	createNQubitOperationQPE(
		numberOfEstimationQubits::Int, 
		numberOfTargetQubits::Int, 
		U::UnitaryOperation, 
		settings::Settings
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Quantum Phase Estimation (QPE) unitary operation** acting on the specified `numberOfEstimationQubits` and `numberOfTargetQubits`, for the given unitary operation ``U`` using qubit indexing `settings`.

The **QPE unitary operation** estimates the eigenphase(s) of the unitary operation ``U``.

Under the hood, it builds the full quantum circuit implementing the **QPE algorithm** and compiles it into a single unitary matrix representing the **QPE unitary operation**.

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (size `numberOfEstimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (size `numberOfTargetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on exactly `numberOfTargetQubits` qubits.
- The resulting unitary operation acts on `numberOfEstimationQubits + numberOfTargetQubits` qubits:
  - The **first `numberOfEstimationQubits`** represent the phase estimation register.
  - The **remaining `numberOfTargetQubits`** represent the target register.
- The qubit indexing is determined by the provided `settings`.

# Arguments
- `numberOfEstimationQubits::Int` — The number of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `numberOfTargetQubits::Int` — The number of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `settings::Settings` — Qubit indexing settings.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **QPE unitary operation**.

# See also
- [`qpeGate!`](@ref)
- [`createNQubitQPEQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)
- [`Settings`](@ref)

# Example
Extract the unitary matrix representation of the **QPE unitary operation** that estimates the eigenphase of a **T-dagger gate** with single-qubit accuracy:
```julia-repl
julia> unitaryOperation = createNQubitOperationQPE(1, 1, createSingleQubitOperationTd(), createSettings(indexType=circuitIndexBigEndian))
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
         1.0-6.12323e-17im       0.0+0.0im       1.01465e-17+6.12323e-17im       0.0+0.0im
         0.0+0.0im          0.853553-0.353553im          0.0+0.0im          0.146447+0.353553im
 1.01465e-17+6.12323e-17im       0.0+0.0im               1.0-1.83697e-16im       0.0+0.0im
         0.0+0.0im          0.146447+0.353553im          0.0+0.0im          0.853553-0.353553im
```
"""
function createNQubitOperationQPE(numberOfEstimationQubits::Int, numberOfTargetQubits::Int, U::UnitaryOperation, settings::Settings)
	quantumCircuit = createNQubitQPEQuantumCircuit(numberOfEstimationQubits, numberOfTargetQubits, U, settings)
	return compileToSingleGate(quantumCircuit)
end

""" 
    createNQubitQPEQuantumCircuit(
		numberOfEstimationQubits::Int, 
		numberOfTargetQubits::Int, 
		U::UnitaryOperation, 
		settings::Settings
	) -> QuantumCircuit{Z, I}

Constructs a **Quantum Phase Estimation (QPE) circuit** that estimates the eigenphase(s) of the unitary operation ``U``.  

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (size `numberOfEstimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (size `numberOfTargetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on exactly `numberOfTargetQubits` qubits.
- The total number of qubits in the circuit is `numberOfEstimationQubits + numberOfTargetQubits`:
  - The **first `numberOfEstimationQubits`** represent the phase estimation register.
  - The **remaining `numberOfTargetQubits`** represent the target register.
- The qubit indexing is determined by the provided `settings`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfEstimationQubits::Int` — The number of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `numberOfTargetQubits::Int` — The number of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `settings::Settings` — Qubit indexing settings.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **QPE**

# See also
- [`createNQubitOperationQPE`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)
- [`Settings`](@ref)

# Example
Create the quantum circuit of a **QPE** that estimates the eigenphase of a **T-dagger gate** with 2-qubit accuracy:
```julia-repl
julia> quantumCircuit = createNQubitQPEQuantumCircuit(2, 1, createSingleQubitOperationTd(), createSettings(indexType=circuitIndexBigEndian))
julia> quantumCircuit.circuit
5-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H"), UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([2], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [2, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [1, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [1, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#57#58"{Vector{Int64}}}(QubiSim.var"#57#58"{Vector{Int64}}([1, 2]), [1, 2], "IQFT")])
```
"""
function createNQubitQPEQuantumCircuit(numberOfEstimationQubits::Int, numberOfTargetQubits::Int, U::UnitaryOperation, settings::Settings)
	numberOfQubits = numberOfEstimationQubits + numberOfTargetQubits
	estimationQubits = Vector(1:numberOfEstimationQubits)
	targetQubits = Vector((numberOfEstimationQubits + 1):numberOfQubits)

	quantumCircuit = createQuantumCircuit(numberOfQubits)
	for m in estimationQubits
        hGate!(quantumCircuit, m)
	end
	estimationQubitsIterator = getEstimationQubitOrder(settings, estimationQubits)

    for m in estimationQubitsIterator
		powerOfU = getPowerOfU(settings, m, numberOfEstimationQubits)
        for _ in 1:powerOfU
            controlledUGate!(quantumCircuit, [m], targetQubits, U)
        end
    end
    iqftGate!(quantumCircuit, getIQFTOrder(settings, estimationQubits))
    return quantumCircuit
end

getEstimationQubitOrder(::Settings{Z, I}, qubits) where {Z, I} = qubits
getEstimationQubitOrder(::Settings{Z, circuitIndexBigEndian}, qubits) where {Z} = reverse(qubits)

getPowerOfU(::Settings{Z, I}, m, numberOfEstimationQubits) where {Z, I} = 2^(m - 1)
getPowerOfU(::Settings{Z, circuitIndexBigEndian}, m, numberOfEstimationQubits) where {Z} = 2^(numberOfEstimationQubits - m)

getIQFTOrder(::Settings{Z, I}, qubits) where {Z, I} = reverse(qubits)
getIQFTOrder(::Settings{Z, circuitIndexBigEndian}, qubits) where {Z} = qubits

#

""" 
	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits; 
		forgetOutcome::Bool = false
	) -> Nothing

	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		sigmas::Sigmas; 
		forgetOutcome::Bool = false
	) -> Nothing

	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		krausOperators::KrausOperators; 
		forgetOutcome::Bool = false
	) -> Nothing

Adds a **measurement gate** to the `quantumCircuit` on the specified `qubits`, with associated measurement direction operators `sigmas` or Kraus operators `krausOperators`. The `forgetOutcome` flag determines whether the measurement outcome is observed or not revealed.

This function performs a quantum measurement based on the provided configuration:
- **Projective measurement (PVM)** is performed when `krausOperators` are not defined.
- **Generalized measurement (POVM)** is performed when `krausOperators` are defined.

Internally, the function constructs a MeasureGate using the given measurement configuration and adds it to the circuit.
During compilation, any **PVM** is automatically converted into a **POVM** by generating the corresponding Kraus operators.

# Arguments
- `quantumCircuit::QuantumCircuit` - The circuit to which the measurement gate is added.
- `qubits::Qubits` - Indices of qubits to be measured, interpreted using the circuit’s `settings`.
- `sigmas::Sigmas` (optional) - Associated measurement direction operators for each qubit, defaulting to ``σ_z``. Used for PVM type measurement.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the generalized POVM type measurement.
- `forgetOutcome::Bool` (optional) — If `true`, the measurement result is not revealed (equivalent to taking the partial trace). If `false`, it is observed. Defaulting to `false`.
  - `false` - Simulates measurement with revealing outcome by collapsing the state to one of the possible outcomes.
  - `true` - Simulates measurement without revealing outcome by creating a mixed state over all possible outcomes.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See Also
- [`QuantumCircuit`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Qubits`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Settings`](@ref)

# Examples

**Example 1: Basic PVM measurement on qubits 1 and 2 along default direction ``σ_z`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[], nothing, false)
```

**Example 2: PVM measurement along custom directions ``σ_x`` and ``σ_y`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2], [sigmaX(), sigmaY()])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[[0.7071067811865476 + 0.0im 0.7071067811865475 - 0.0im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im], [0.7071067811865476 + 0.0im 4.329780281177466e-17 - 0.7071067811865475im; 4.329780281177466e-17 + 0.7071067811865475im -0.7071067811865476 + 0.0im]], nothing, false)
```

**Example 3: PVM measurement along direction ``σ_z`` with outcome NOT revealed to the observer**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2], [sigmaZ(), sigmaZ()], forgetOutcome = true)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 0.0im], [1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 0.0im]], nothing, true)
```

**Example 4: POVM measurement using Kraus operators for a PVM measurement along direction ``σ_z`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> measureGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]), forgetOutcome = false)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.POVM, [1], nothing, KrausOperator[KrausOperator(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], "0"), KrausOperator(ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], "1")], false)
```
"""
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), Sigmas(), forgetOutcome))
end
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, sigmas::Sigmas; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), sigmas, forgetOutcome))
end
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, krausOperators::KrausOperators; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), krausOperators, forgetOutcome))
end

#

""" 
	quantumChannelGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		krausOperators::KrausOperators
	) -> Nothing

Adds a **quantum channel gate** (also known as a *completely positive trace-preserving map*, or CPTP map) to the `quantumCircuit` on the specified `qubits` using the provided `krausOperators`.

This function models a **superoperator** — a general type of quantum transformation that includes noise, decoherence, or collapse due to measurement.

Internally, it creates a `QuantumChannelGate` with the given Kraus operators and adds it to the circuit.

# Arguments
- `quantumCircuit::QuantumCircuit` - The circuit to which the quantum channel gate is added.
- `qubits::Qubits` - Indices of the qubits the gate acts on, interpreted using the circuit’s `settings`.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the quantum channel.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See Also
- [`QuantumCircuit`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`Qubits`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Settings`](@ref)

# Example
Add a **quantum channel gate** to a single-qubit circuit using Kraus operators for a projective measurement along direction ``σ_z`` with the outcome not revealed:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
QuantumChannelGate([1], Matrix{ComplexF64}[[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]])
```
"""
function quantumChannelGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, krausOperators::KrausOperators)
	barrier!(quantumCircuit)
	addGate!(quantumCircuit, QuantumChannelGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), krausOperators))
	barrier!(quantumCircuit)
end

#

""" 
    createSettings(;
		zeroBasedNumbering = false, 
		indexType = circuitIndexBigEndian
	) -> Settings{Z, I}

Creates a compile-time structure defining **global indexing settings** for qubits in a quantum circuit.  
This type encodes its parameters as type-level values for efficient specialization and dispatch.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `zeroBasedNumbering::Bool`: 
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1.

- `indexType::IndexType`: 
  - Determines how qubits are mapped between circuit view and byte representation. 
  - Can be `circuitIndexBigEndian`, `circuitIndexLittleEndian`, or `byteIndex`.

# Returns
- `Settings{Z, I}` — Creates a specialized instance, where both parameters `Z` and `I` are encoded at compile time.

# See Also
- [`IndexType`](@ref)
- [`Settings`](@ref)

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# Example
Create the default indexing settings (zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
```julia-repl
julia> settings = createSettings()
Settings{false, circuitIndexBigEndian}()
```
"""
function createSettings(;zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
	Settings(zeroBasedNumbering, indexType)
end

""" 
    createQuantumCircuitWithSettings(
		numberOfQubits::Int; 
		circuit = Circuit(), 
		settings = createSettings()
	) -> QuantumCircuit{Z, I}


Creates a `QuantumCircuit` with the specified `numberOfQubits`, an optional initial `circuit`, using custom qubit indexing `settings`.

This function initializes a quantum circuit structure. By default, it uses an empty circuit and default settings.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` - The total number of qubits in the circuit.
- `circuit::Circuit` (optional) - A sequence of circuit steps, each being a subtype of `Step` (`UnitaryStep`, `MeasurementStep`, `QuantumChannelStep`). Defaults to an empty `Circuit()`.
- `settings::Settings` (optional) - Qubit indexing settings. Defaults to `createSettings()`.

# Returns
- `QuantumCircuit{Z, I}` - A quantum circuit that includes the provided number of qubits, initialization from the given `circuit`, using indexing settings.

# See Also
- [`QuantumCircuit`](@ref)
- [`Circuit`](@ref)
- [`Settings`](@ref)
- [`createSettings`](@ref)

# Example
Create an empty quantum circuit with 2 qubits using default settings:
```julia-repl
julia> qc = createQuantumCircuitWithSettings(2)
QuantumCircuit{false, circuitIndexBigEndian}(Settings{false, circuitIndexBigEndian}(), 2, [1, 1], Step[])
```
"""
function createQuantumCircuitWithSettings(numberOfQubits::Int; circuit = Circuit(), settings = createSettings())
	# set numberOfStepsOnQubits to length(circuit)==number of steps in circuit as barrier call
	QuantumCircuit(settings, numberOfQubits, maximum([1,length(circuit)])*ones(Int, numberOfQubits), circuit)
end

""" 
    createQuantumCircuit(
		numberOfQubits::Int; 
		circuit = Circuit(), 
		zeroBasedNumbering = false, 
		indexType = circuitIndexBigEndian
	) -> QuantumCircuit{Z, I}

Creates a `QuantumCircuit` with the specified `numberOfQubits`, an optional initial `circuit`, using indexing settings `zeroBasedNumbering` and `indexType`.

This function initializes a quantum circuit structure. By default, it uses an empty circuit and default settings.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` - The total number of qubits in the circuit.
- `circuit::Circuit` (optional) - A sequence of circuit steps, each being a subtype of `Step` (`UnitaryStep`, `MeasurementStep`, `QuantumChannelStep`). Defaults to an empty `Circuit()`.
- `zeroBasedNumbering::Bool` (optional): 
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1. Defaults to `false`.

- `indexType::IndexType` (optional): 
  - Determines how qubits are mapped between circuit view and byte representation. 
  - Can be `circuitIndexBigEndian`, `circuitIndexLittleEndian`, or `byteIndex`. Defaults to `circuitIndexBigEndian`.

# Returns
- `QuantumCircuit{Z, I}` - A quantum circuit that includes the provided number of qubits, initialization from the given `circuit`, using indexing settings.

# See Also
- [`QuantumCircuit`](@ref)
- [`Circuit`](@ref)
- [`IndexType`](@ref)
- [`Settings`](@ref)

# Example
Create an empty quantum circuit with 2 qubits using default settings:
```julia-repl
julia> qc = createQuantumCircuit(2)
QuantumCircuit{false, circuitIndexBigEndian}(Settings{false, circuitIndexBigEndian}(), 2, [1, 1], Step[])
```
"""
function createQuantumCircuit(numberOfQubits::Int; circuit = Circuit(), zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
	createQuantumCircuitWithSettings(numberOfQubits; circuit = circuit, settings = Settings(zeroBasedNumbering, indexType))
end

""" 
    concatenateQuantumCircuits(
		a::QuantumCircuit, 
		b::QuantumCircuit
	) -> QuantumCircuit{Z, I}

Concatenates two quantum circuits `a` and `b` into a single `QuantumCircuit`.

Both circuits must have the same number of qubits and identical settings. The resulting circuit includes all operations from `a` followed by those from `b`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `a::QuantumCircuit` - The first quantum circuit.
- `b::QuantumCircuit` - The second quantum circuit to append after `a`.

# Returns
- `QuantumCircuit{Z, I}` - A new quantum circuit that combines the operations from `a` and `b`.

# Throws
- A `String` exception if `a` and `b` have different `numberOfQubits` or mismatched settings.

# See Also
- [`QuantumCircuit`](@ref)
- [`Settings`](@ref)

# Example
Concatenate a circuit with a **Hadamard gate** and one with an **Pauli-X gate**:
```julia-repl
julia> qca = createQuantumCircuit(1)
julia> hGate!(qca, 1)
julia> qcb = createQuantumCircuit(1)
julia> xGate!(qcb, 1)
julia> qc = concatenateQuantumCircuits(qca, qcb)
QuantumCircuit{false, circuitIndexBigEndian}(Settings{false, circuitIndexBigEndian}(), 1, [3], Step[UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")]), UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [1], "X")])])
```
"""
function concatenateQuantumCircuits(a::QuantumCircuit, b::QuantumCircuit)
	if a.numberOfQubits == b.numberOfQubits && a.settings == b.settings
		c = createQuantumCircuitWithSettings(a.numberOfQubits, settings = a.settings)
		for gate in a.circuit
			push!(c.circuit, gate)
		end
		for gate in b.circuit
			push!(c.circuit, gate)
		end
		# set numberOfStepsOnQubits to length(circuit)==number of steps in circuit as barrier call
		c.numberOfStepsOnQubits = (maximum([1,length(c.circuit)])+1)*ones(Int64, c.numberOfQubits)
		return c
	else
		throw("Cannot concatenate quantum circuits: numberOfQubits or settings do not match.");
	end
end

""" 
	inverseQuantumCircuit(
		quantumCircuit::QuantumCircuit{Z, I}
	) -> QuantumCircuit{Z, I}

Inverses the given `quantumCircuit`.

This function constructs a new `QuantumCircuit` by reversing the order of steps in the input circuit. 

It assumes that the gates used are unitary and that applying them in reverse order effectively inverts the computation.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `quantumCircuit::QuantumCircuit{Z, I}` - The input quantum circuit to be inverted.

# Returns
- `QuantumCircuit{Z, I}` - A new quantum circuit that represents the inverse of the input circuit.

# See Also
- [`QuantumCircuit`](@ref)
- [`Circuit`](@ref)

# Example
Invert a quantum circuit that applies a **Hadamard gate** followed by an **Pauli-X gate**:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> inverseQuantumCircuit(qc)
QuantumCircuit{false, circuitIndexBigEndian}(Settings{false, circuitIndexBigEndian}(), 1, [3], Step[UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [1], "X")]), UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])])
```
"""
function inverseQuantumCircuit(quantumCircuit::QuantumCircuit)
    numberOfQubits = quantumCircuit.numberOfQubits
    settings = quantumCircuit.settings
    steps = quantumCircuit.circuit

    inversedQuantumCircuit = createQuantumCircuitWithSettings(numberOfQubits, settings = settings)
    for step in reverse(steps)
        push!(inversedQuantumCircuit.circuit, step)
    end

    inversedQuantumCircuit.numberOfStepsOnQubits = (length(steps) + 1) * ones(Int64, numberOfQubits)
    return inversedQuantumCircuit
end

""" 
    compileToSingleGate(
		quantumCircuit::QuantumCircuit
	) -> UnitaryOperation

Compiles a given `quantumCircuit` that consist of only `UnitaryGate` types into a single `UnitaryOperation`.

# Arguments
- `quantumCircuit::QuantumCircuit` - The quantum circuit to compile.

# Returns
- `UnitaryOperation` - A single unitary representation of the quantum circuit.

# See Also
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`UnitaryOperation`](@ref)

# Throws
- A `String` exception if the quantum circuit does not contain only `UnitaryGate` types.

# Example

Step 1 - Compile a quantum circuit that applies a **Hadamard gate** followed by an **Pauli-X gate** into a single **HX unitary operation**:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> hxGate = compileToSingleGate(qc)
julia> hxGate.U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im  -0.707107+0.0im
 0.707107+0.0im   0.707107+0.0im
```
Step 2 - Create and compile a quantum circuit containing the single **HX unitary operation* created above:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> unitaryUGate!(qc, [1], hxGate)
julia> qp = compileQuantumCircuit(qc)
julia> qp.program[1].U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im  -0.707107+0.0im
 0.707107+0.0im   0.707107+0.0im
```
"""
function compileToSingleGate(quantumCircuit::QuantumCircuit)
	quantumProgram = compileQuantumCircuit(quantumCircuit; optimizeNumberOfSteps = true)
	if length(quantumProgram.program) == 1
		return quantumProgram.program[1]
	else
		throw("The quantumCircuit does not have only UnitaryGate types")
	end
end

""" 
    tensorProduct(
		a, 
		b
	) -> AbstractArray

Computes the Kronecker (tensor) product ``a ⊗ b`` of two quantum operations or states ``a`` and ``b``.

This operation is used to combine two quantum systems into a single composite system. It supports both vectors (quantum states) and matrices (quantum operations).

# Arguments
- `a` - A vector or matrix ``a`` representing a quantum state or operation.
- `b` - A vector or matrix ``b`` representing a quantum state or operation.

# Returns
- `AbstractArray` - A new array representing the tensor (Kronecker) product of ``a`` and ``b``. 

# Examples

**Example 1: Tensor product of the ``|0⟩`` and ``|+⟩`` vector states**
```julia-repl
julia> tensorProduct([1; 0], [1; 1]/sqrt(2))
4-element Vector{Float64}:
 0.7071067811865475
 0.7071067811865475
 0.0
 0.0
```
**Example 2: Tensor product of the identity and Hadamard unitary operations**
```julia-repl
julia> tensorProduct([1 0; 0 1], [1 1; 1 -1]/sqrt(2))
4×4 Matrix{Float64}:
 0.707107   0.707107  0.0        0.0
 0.707107  -0.707107  0.0       -0.0
 0.0        0.0       0.707107   0.707107
 0.0       -0.0       0.707107  -0.707107
```
**Example 3: Tensor product of the Hadamard and identity unitary operations**
```julia-repl
julia> tensorProduct([1 1; 1 -1]/sqrt(2), [1 0; 0 1])
4×4 Matrix{Float64}:
 0.707107  0.0        0.707107   0.0
 0.0       0.707107   0.0        0.707107
 0.707107  0.0       -0.707107  -0.0
 0.0       0.707107  -0.0       -0.707107
```
"""
function tensorProduct(a, b)
	kron(a, b)
end

"""
    createToggleSwapList(
		numberOfQubits::Int
	) -> Vector{Int}

Generates a list of indices representing bit-reversed ordering for a quantum system with `numberOfQubits` qubits.

# Arguments
- `numberOfQubits::Int` - The number of qubits ``M``.

# Returns
- `Vector{Int}` - A list of indices (1-based) of length ``2^M`` that perform bit-reversal reordering.

# Notes
- The returned list can be used to permute state vectors or matrices: `UnitaryOperation.U[swapList, swapList]`, `VectorState.q[swapList]` and `DensityState.rho[swapList, swapList]`.

# See Also
- [`UnitaryOperation`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1: Reordering for 2 qubits**

Binary representations:
- index `1` → `0` (00) → (00) `0` → index `1`
- index `2` → `1` (01) → (10) `2` → index `3`
- index `3` → `2` (10) → (01) `1` → index `2`
- index `4` → `3` (11) → (11) `3` → index `4`

```julia-repl
julia> createToggleSwapList(2)
4-element Vector{Int64}:
 1
 3
 2
 4
```

**Example 2: Reordering for 3 qubits**

Binary representations:
- index `1` → `0` (000) → (000) `0` → index `1`
- index `2` → `1` (001) → (100) `4` → index `5`
- index `3` → `2` (010) → (010) `2` → index `3`
- index `4` → `3` (011) → (110) `6` → index `7`
- index `5` → `4` (100) → (001) `1` → index `2`
- index `6` → `5` (101) → (101) `5` → index `6`
- index `7` → `6` (110) → (011) `3` → index `4`
- index `8` → `7` (111) → (111) `7` → index `8`

```julia-repl
julia> createToggleSwapList(3)
8-element Vector{Int64}:
 1
 5
 3
 7
 2
 6
 4
 8
```
"""
function createToggleSwapList(numberOfQubits::Int)
	n = 2^numberOfQubits
	swapList = ones(Int, n)
	for m in 1:n
        # Compute bit-reversed index (0-based), convert to 1-based
  		# bits = reverse(digits(m - 1, base = 2, pad = numberOfQubits)')
  		# swapList[m] = sum([bits[i] * 2^(i - 1) for i in 1:numberOfQubits]) + 1
  		bits = reverse(digits(m - 1, base = 2, pad = numberOfQubits))
  		swapList[m] = sum(bits[i] * 2^(i - 1) for i in 1:numberOfQubits) + 1
	end
	return swapList
end

function constructKrausOperatorsForPVMMeasurement(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas)
	list = generateMeasurementIndexMap(numberOfQubits, qubits)
	transform = buildTransformMatrix(numberOfQubits, qubits, sigmas)
	labels = generateMeasurementLabels(numberOfQubits, qubits)

	# Sort list/labels for little-endian order
	ind = sortperm(list[:, 1])
	list = list[ind, :]
	labels = labels[ind]

	return constructKrausOperators(numberOfQubits, list, transform, labels)
end

function generateMeasurementIndexMap(numberOfQubits::Int, qubits::Qubits)
	numMeasured = length(qubits)
	numOutcomes = 2^numMeasured
	list = zeros(Int, numOutcomes, 2^(numberOfQubits - numMeasured))

	for m in 0:numOutcomes-1
		targetBits = reverse(digits(m, base = 2, pad = numMeasured)')
		k = 1
		for n in 0:(2^numberOfQubits - 1)
			fullBits = reverse(digits(n, base = 2, pad = numberOfQubits)')
			if isequal(fullBits[qubits]', targetBits)
				list[m+1, k] = n + 1
				k += 1
			end
		end
	end

	return list
end

function buildTransformMatrix(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas)
	transform = 1
	for j in 1:numberOfQubits
		ind = indexin(j, qubits)[1]
		V = isnothing(ind) ? Matrix{Float64}(I, 2, 2) : sigmas[ind]
		transform = tensorProduct(transform, V)
	end
	return transform
end

function generateMeasurementLabels(numberOfQubits::Int, qubits::Qubits)
	numMeasured = length(qubits)
	numOutcomes = 2^numMeasured
	labels = String[]

	for r in 1:numOutcomes
		label = []
		for k in 1:numberOfQubits
			push!(label, '*')
		end
		s = reverse(digits(r - 1, base = 2, pad = numMeasured)')
		for (k, qubit) in enumerate(qubits)
			label[qubit] = s[k]
		end
		push!(labels, join(string.(label)))
	end

	return labels
end

function constructKrausOperators(numberOfQubits::Int, list::Matrix{Int}, transform::AbstractMatrix, labels::Vector{String})
	N = 2^numberOfQubits
	krausOperators = KrausOperators()

	for m in axes(list, 1)
		Mm = zeros(N, N)
		for n in axes(list, 2)
			Mm[list[m, n], list[m, n]] = 1
		end
		Mm = transform' * Mm * transform # Mm = T'*|m><m|*T
		Em = Mm # Em = T'*|m><m|*T
		push!(krausOperators, KrausOperator(Em, labels[m]))
	end

	return krausOperators
end

""" 
    sigmaX() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the X-axis on the Bloch sphere.

This operator defines a projective measurement along the X direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaX()`, the measurement direction corresponds to spherical polar angles ``θ = π/2`` and ``φ = 0``, aligning with the positive X-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the X-axis.

# Calling structure
```julia-repl
sigmaX()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
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

# Example
```julia-repl
julia> sigmaX()
2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107-0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
function sigmaX()
	return sigmaN(pi/2, 0.)
	# sx=pauliX() # [0 1;1 0]
	# V=eigvecs(sx)
	# D=eigvals(sx)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# V[:,2]=-V[:,2]
	# return V
end

""" 
    sigmaY() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the Y-axis on the Bloch sphere.

This operator defines a projective measurement along the Y direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaY()`, the measurement direction corresponds to spherical polar angles ``θ = π/2`` and ``φ = π/2``, aligning with the positive Z-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the Y-axis.

# Calling structure
```julia-repl
sigmaY()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
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

# Example
```julia-repl
julia> sigmaY()
2×2 Matrix{ComplexF64}:
    0.707107+0.0im   0.0-0.707107im
    0.0+0.707107im  -0.707107+0.0im
```
"""
function sigmaY()
	return sigmaN(pi/2, pi/2)
	# sy=pauliY() # [0 -1im;1im 0]
	# V=eigvecs(sy)
	# D=eigvals(sy)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# V=1im*V
	# return V
end

""" 
    sigmaZ() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the Z-axis on the Bloch sphere.

This operator defines a projective measurement along the Z direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaZ()`, the measurement direction corresponds to spherical polar angles ``θ = 0`` and ``φ = 0``, aligning with the positive Z-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the Z-axis.

# Calling structure
```julia-repl
sigmaZ()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaY`](@ref)
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

# Example
```julia-repl
julia> sigmaZ()
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0-0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function sigmaZ()
	return sigmaN(0., 0.)
	# sz=pauliZ() # [1 0;0 -1]
	# V=eigvecs(sz)
	# D=eigvals(sz)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# return V
	# return sz
end

""" 
    sigmaN(
		n::Vector{Float64}
	) -> ComplexMatrix
    
	sigmaN(
		theta::Float64, 
		phi::Float64
	) -> ComplexMatrix

Creates the **measurement direction operator** corresponding to an arbitrary direction on the Bloch sphere, specified either by a 3D unit vector ``n`` or by spherical polar angles ``θ`` and ``φ``.

Each **measurement direction operator** can be parameterized by a 3D unit vector ``n = [n_x, n_y, n_z]`` on the Bloch sphere:
```math
    M(n) = n ⋅ σ = n_x ⋅ σ_x + n_y ⋅ σ_y + n_z ⋅ σ_z
```

with ``σ = [σ_x, σ_y, σ_z]`` the vector of Pauli matrices associated with the X, Y, and Z axes, respectively.

Alternatively, if the unit vector ``n`` is expressed in terms of spherical polar angles ``θ`` and ``φ``, the operator takes the form:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

# Arguments
- `n::Vector{Float64}` - A 3-dimensional unit vector ``n = [n_x, n_y, n_z]`` specifying a direction on the Bloch sphere.
- `theta::Float64` - Polar angle ``θ`` in radians (angle from the Z-axis).
- `phi::Float64` - Azimuthal angle ``φ`` in radians (angle from the X-axis in the XY-plane).

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** for the specified direction.

# Calling structure
```julia-repl
sigmaN(theta::Float64, phi::Float64)
    └── sigmaN(n::Vector{Float64})
        └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaY`](@ref)
- [`sigmaZ`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

# Examples

**Example 1: Using a unit vector direction**

```julia-repl
julia> sigmaN(vec([1 1 1])/sqrt(3))
2×2 Matrix{ComplexF64}:
 0.57735+0.0im       0.57735-0.57735im
 0.57735+0.57735im  -0.57735+0.0im
```

**Example 2: Using spherical coordinates ``θ = π/3``, ``φ = π/4``**

```julia-repl
julia> sigmaN(pi/3., pi/4.)
2×2 Matrix{ComplexF64}:
 0.866025+0.0im        0.353553-0.353553im
 0.353553+0.353553im  -0.866025+0.0im
```
"""
function sigmaN(theta::Float64, phi::Float64)
	return sigmaN(vec([sin(theta/2)*cos(phi) sin(theta/2)*sin(phi) cos(theta/2)]))
end

function sigmaN(n::AbstractVector{<:Number})
	if length(n) == 3
		return createSingleQubitMeasureOperator(Float64.(n))
	else
		throw("length of vector n is not equal to three")
	end
	# V=eigvecs(sn)
	# D=eigvals(sn)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# return V
end

# # Calling structure
# ```julia-repl
# createSingleQubitMeasureOperator(theta::Float64, phi::Float64)
#     └── createSingleQubitMeasureOperator(n::Vector{Float64})
#         └── [matrix construction]
# ```

function createSingleQubitMeasureOperator(theta::Float64, phi::Float64)
	return createSingleQubitMeasureOperator(vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]))
end

function createSingleQubitMeasureOperator(n::Vector{Float64})
	# nsigma = a*sigma with sigma the 3 pauli matrices
	# rho = n(1)*[[0,1];[1,0]]+n(2)*[[0,-i];[i,0]]+n(3)*[[1,0];[0,-1]]
	return [n[3] n[1]-im*n[2]; n[1]+im*n[2] -n[3]]
end

""" 
    pauliX() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-X unitary operation**.

The **Pauli-X unitary operation** flips the bit of a single qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

# Example
```julia-repl
julia> pauliX()
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
function pauliX()
	return ComplexMatrix([0 1;1 0])
end

""" 
    pauliY() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-Y unitary operation**.

The **Pauli-Y unitary operation** flips both the bit and the phase of a single qubit.

# Example
```julia-repl
julia> pauliY()
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
function pauliY()
	return ComplexMatrix([0 -1im;1im 0])
end

""" 
    pauliZ() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-Z unitary operation**.

The **Pauli-Z unitary operation** flips the phase of a single qubit.

# Example
```julia-repl
julia> pauliZ()
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function pauliZ()
	return ComplexMatrix([1 0;0 -1])
end

""" 
    createSingleQubitBlochDensityState(
		n::Vector{Float64}
	) -> DensityState
	
	createSingleQubitBlochDensityState(
		theta::Float64, 
		phi::Float64
	) -> DensityState

Creates a single qubit **Bloch density state** in the direction specified either by a 3D unit vector ``n`` or by spherical polar angles ``θ`` and ``φ``.

The Bloch density state is given by: 
```math
ρ = 0.5 * (I + n ⋅ σ) = 0.5 * (I + n_x ⋅ σ_x + n_y ⋅ σ_y + n_z ⋅ σ_z)
```

where:
- ``I`` is the 2×2 identity matrix,
- ``n = [n_x, n_y, n_z]`` is a real 3D unit vector representing a point on the Bloch sphere,
- ``n = [sin(θ)⋅cos(φ), sin(θ)⋅sin(φ), cos(θ)]`` when parameterized by spherical polar angles,
- ``σ = [σ_x, σ_y, σ_z]`` is the vector of Pauli matrices associated with the X, Y, and Z axes, respectively.

# Arguments
- `n::Vector{Float64}` - A 3-dimensional unit vector ``n = [n_x, n_y, n_z]`` specifying a direction on the Bloch sphere.
- `theta::Float64` - Polar angle ``θ`` in radians (angle from the Z-axis).
- `phi::Float64` - Azimuthal angle ``φ`` in radians (angle from the X-axis in the XY-plane).

# Returns
- `DensityState` - A single-qubit **Bloch density state** ρ aligned with the specified Bloch vector direction.

# Notes
- The resulting matrix ``ρ`` is Hermitian, positive semi-definite, and has trace 1.
- This representation is useful for simulating qubit states in both pure and mixed form.
- When ``‖n‖ = 1``, the state is pure. When ``‖n‖ < 1``, the state is mixed.

# Calling structure
```julia-repl
createSingleQubitBlochDensityState(theta::Float64, phi::Float64)
    └── createSingleQubitBlochDensityState(n::Vector{Float64})
        └── DensityState(...)
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1: Using a normalized Bloch vector**

```julia-repl
julia> createSingleQubitBlochDensityState(vec([1 1 1])/sqrt(3)).rho
2×2 Matrix{ComplexF64}:
 0.788675+0.0im       0.288675-0.288675im
 0.288675+0.288675im  0.211325+0.0im
```

**Example 2: Using spherical coordinates ``θ = π/3``, ``φ = π/4``**

```julia-repl
julia> createSingleQubitBlochDensityState(pi/3., pi/4.).rho
2×2 Matrix{ComplexF64}:
     0.75+0.0im       0.306186-0.306186im
 0.306186+0.306186im      0.25+0.0im
```
"""
function createSingleQubitBlochDensityState(n::Vector{Float64})
	# rho = 0.5*(identity+n*sigma) with sigma the 3 pauli matrices
	# rho = 0.5*(eye(2)+n(1)*[[0,1];[1,0]]+n(2)*[[0,-i];[i,0]]+n(3)*[[1,0];[0,-1]]);
	return DensityState(0.5*[1+n[3] n[1]-im*n[2]; n[1]+im*n[2] 1-n[3]])
end

function createSingleQubitBlochDensityState(theta::Float64, phi::Float64)
	return createSingleQubitBlochDensityState(vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]))
end

function cascadeOperations(U1::UnitaryOperation, U2::UnitaryOperation)
	UnitaryOperation(U2.U * U1.U) # watch out for reverse order!!! first U1 and then U2, i.e., U2(U1(q))=U2*U1*q
end

function mergeTwoOperators(A::ComplexMatrix, B::ComplexMatrix, listA::Vector{Int64}, listB::Vector{Int64})
	nrBitsA = log2(size(A, 1))
	nrBitsB = log2(size(B, 1))
	nrBitsU = Int(nrBitsA + nrBitsB)
	N = 2^nrBitsU
	U = zeros(N, N)
	iA = ones(Int, N) # to force a Vector that can be used to index into Matrix
	if !isempty(listA)
		for k in 0:(N-1)
			for (m, element) in enumerate(listA)
				iA[k+1] += mod(floor(k / 2^(nrBitsU - element)), 2) * 2^(length(listA) - m);
			end
		end
	end
	iB = ones(Int, N) # to force a Vector that can be used to index into Matrix
	if !isempty(listB)
		for k in 0:(N-1)
			for (m, element) in enumerate(listB)
				iB[k+1] += mod(floor(k / 2^(nrBitsU - element)), 2) * 2^(length(listB) - m);
			end
		end
	end
	U = A[iA, iA] .* B[iB, iB]
end

function moveNQubitOperationToListOfQubits(U::UnitaryOperation, list::Vector{Int64}, n::Int)
	if size(U.U,2) == size(U.U,1) == 2^length(list)
		U = UnitaryOperation(mergeTwoOperators(U.U, createNQubitOperatorId(n-length(list)), list, setdiff(1:n, list)))
	else
		throw("size of UnitaryOperation does not match length of qubits list")
	end
end

function moveNQubitKrausOperatorToListOfQubits(E::KrausOperator, list::Vector{Int64}, n::Int)
	if size(E.E,2) == size(E.E,1) == 2^length(list)
		E = KrausOperator(mergeTwoOperators(E.E, createNQubitOperatorId(n-length(list)), list, setdiff(1:n, list)), E.label)
	else
		throw("size of KrausOperator does not match length of qubits list")
	end
end

@enum StepType unitary measure quantumChannel unknown

""" 
    compileQuantumCircuit(
		quantumCircuit::QuantumCircuit; 
		optimizeNumberOfSteps::Bool = false
	) -> QuantumProgram{Z, I}

Compiles a `quantumCircuit` into a `QuantumProgram` by translating each logical step into a corresponding quantum operation.

Each step in `quantumCircuit.circuit::Vector{Step}` is a vector of gates, where **all gates must be of the same gate type**. Supported gate types include:
- `UnitaryGate`
- `MeasureGate`
- `QuantumChannelGate`

### Step Compilation Behavior

- **Unitary Steps**:  
  For steps containing only `UnitaryGate`s, a tensor product of all gate unitary matrices in that step is constructed to form a single `UnitaryOperation`.

- **Measurement Steps**:  
  For steps containing `MeasureGate`s, a `MeasureOperation` or `MeasureAndForgetOperation` is created based on the qubits, measurement direction operators (`sigmas`) or Kraus operators (`krausOperators`), and measurement settings.

- **Quantum Channel Steps**:  
  For steps containing `QuantumChannelGate`s, the defined Kraus operators are used to form a `QuantumChannelOperation`.

### Optimization

If `optimizeNumberOfSteps = true`, all **consecutive steps** made up entirely of `UnitaryGate`s are merged into a **single** `UnitaryOperation` to reduce circuit depth and improve execution performance.

Otherwise (default behavior), each step is compiled into a **separate** `UnitaryOperation`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `quantumCircuit::QuantumCircuit` - The quantum circuit to be compiled.
- `optimizeNumberOfSteps::Bool` (optional) - Whether to combine consecutive unitary steps into one operation. Defaults to `false`.

# Returns
- `QuantumProgram{Z, I}`: The compiled quantum program ready for execution.

# See Also
- [`QuantumCircuit`](@ref)
- [`QuantumProgram`](@ref)
- [`Step`](@ref)
- [`UnitaryGate`](@ref)
- [`MeasureGate`](@ref)
- [`QuantumChannelGate`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)

# Examples

**Example 1: Compile a quantum circuit with separate Hadamard and Pauli-X gates**

```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> xGate!(qc, 1)
julia> qc.circuit
2-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [1], "X")])
julia> qp = compileQuantumCircuit(qc)
julia> qp.program
2-element Vector{Operation}:
 UnitaryOperation(ComplexF64[0.7071067811865476 + 0.0im 0.7071067811865475 - 8.659560562354932e-17im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 8.659560562354934e-17im])
 UnitaryOperation(ComplexF64[6.123233995736766e-17 + 0.0im 1.0 - 1.2246467991473532e-16im; 1.0 + 0.0im -6.123233995736766e-17 + 7.498798913309288e-33im])
```

**Example 2: Compile the same circuit with step optimization**

```julia-repl
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> qp.program
1-element Vector{Operation}:
 UnitaryOperation(ComplexF64[0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im; 0.7071067811865476 + 0.0im 0.7071067811865475 + 0.0im])
```
"""
function compileQuantumCircuit(quantumCircuit::QuantumCircuit; optimizeNumberOfSteps::Bool = false)
	quantumProgram = QuantumProgram(quantumCircuit.settings, quantumCircuit.numberOfQubits, Program())
    for step in quantumCircuit.circuit
        compileStep!(quantumProgram, step, quantumCircuit.numberOfQubits, optimizeNumberOfSteps)
    end
    return quantumProgram
end

function compileStep!(program::QuantumProgram, step::UnitaryStep, numberOfQubits::Int, optimize::Bool)
    UStep = createNQubitOperationId(numberOfQubits)
    for gate in step.gates
		UGate = gate.unitaryOperationFactory()
        UGate = moveNQubitOperationToListOfQubits(UGate, gate.qubits, numberOfQubits)
        UStep = cascadeOperations(UStep, UGate)
    end

	# addStepToProgram!(program, UStep, Val(optimize))
    if optimize && !isempty(program.program) && program.program[end] isa UnitaryOperation
        # Merge with previous unitary
        program.program[end] = cascadeOperations(program.program[end], UStep)
    else
        push!(program.program, UStep)
    end
end

# function addStepToProgram!(program::QuantumProgram, UStep::UnitaryOperation, ::Val{true})
# 	if !isempty(program.program) && program.program[end] isa UnitaryOperation
#         # Merge with previous unitary
#         program.program[end] = cascadeOperations(program.program[end], UStep)
# 	else
# 	    push!(program.program, UStep)
#     end
# end

# function addStepToProgram!(program::QuantumProgram, UStep::UnitaryOperation, ::Val{false})
#     push!(program.program, UStep)
# end

function compileStep!(program::QuantumProgram, step::MeasurementStep, numberOfQubits::Int, optimize::Bool)
    gate = step.gate
    qubits = gate.qubits
	if gate.measureGateType == PVM
		sigmas = isempty(gate.sigmas) ? [sigmaZ() for _ in qubits] : gate.sigmas
    	krausOperators = constructKrausOperatorsForPVMMeasurement(numberOfQubits, qubits, sigmas)
	elseif gate.measureGateType == POVM
		krausOperators = [moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits) for krausOperator in gate.krausOperators]
	end

	measurementOperators = [krausOperators[m].E'*krausOperators[m].E for m in 1:length(krausOperators)]
	labels=[krausOperator.label for krausOperator in krausOperators]

	if gate.forgetOutcome
	    push!(program.program, MeasureAndForgetOperation(krausOperators, measurementOperators, labels))
	else
    	push!(program.program, MeasureOperation(krausOperators, measurementOperators, labels))
	end
end

function compileStep!(program::QuantumProgram, step::QuantumChannelStep, numberOfQubits::Int, optimize::Bool)
    gate = step.gate
    qubits = gate.qubits
	krausOperators = [moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits) for krausOperator in gate.krausOperators]
    push!(program.program, QuantumChannelOperation(krausOperators))
end

function compileStep!(program::QuantumProgram, step::Step, numberOfQubits::Int, optimize::Bool)
    error("Unsupported quantum step type: $(typeof(step))")
end

"""
    @enum StateSpace vector density

Supported quantum mechanical state representations:

- `vector`: **State vector** representation ``|Ψ⟩``
  Represents pure quantum states as complex-valued vectors in Hilbert space.

- `density`: **Density matrix** representation ``ρ``
  Represents statistical mixtures of pure states and supports both pure and mixed quantum states via a positive semi-definite matrix with unit trace.
"""
@enum StateSpace vector density

"""
    createByteIndexVector(
		value::Int, 
		numberOfQubits::Int
	) -> Matrix{Float64}

Constructs a **multi-qubit product vector state** from an integer `value` for a system of `numberOfQubits` qubits.

The resulting quantum vector state has the form: ``|ψ⟩ = |ψ₁⟩ ⊗ |ψ₂⟩ ⊗ ⋯ ⊗ |ψₙ⟩`` where each individual qubit is in the basis vector state: ``|ψⱼ⟩ = αⱼ|0⟩ + βⱼ|1⟩``

The function returns a 2×`numberOfQubits` matrix of the form: ``[α_1 ... α_n; β_1 ... β_n]``.

Each column represents one qubit, and is initialized to:
- ``|0⟩ = [1; 0]`` if the corresponding bit in the binary representation of `value` is 0
- ``|1⟩ = [0; 1]`` if the bit is 1

The least significant bit of `value` corresponds to the last qubit (rightmost column), consistent with standard quantum computing conventions.

This representation is useful as input for `createInitialQubitState` when using `indexType = byteIndex`.

# Arguments
- `value::Integer` - Integer value to convert into binary qubit states.
- `numberOfQubits::Integer` - Number of qubits to encode.

# Returns
- `Matrix{Float64}` - A 2×`numberOfQubits` matrix representing the qubit product state.

# See Also
- [`IndexType`](@ref)
- [`createInitialQubitState`](@ref)

# Example
Create a 3-qubit state from the value 6 (binary 110), producing ``|ψ⟩ = |1⟩ ⊗ |1⟩ ⊗ |0⟩``:
```julia-repl
julia> createByteIndexVector(6, 3)
2×3 Matrix{Float64}:
 0.0  0.0  1.0
 1.0  1.0  0.0
```
"""
function createByteIndexVector(value::Int, numberOfQubits::Int)
    valueBinary = reverse(digits(value, base = 2, pad = numberOfQubits)')
    byteIndexVector = zeros(2, numberOfQubits)
    for k in 1:numberOfQubits
        if valueBinary[k] == 0
            byteIndexVector[1, k] = 1.
        else
            byteIndexVector[2, k] = 1.
        end
    end
    return byteIndexVector
end

"""
	createInitialQubitState(
		stateSpace::StateSpace, 
		value; 
		indexType = circuitIndexBigEndian, 
		blochRepresentation = false
	) -> State

Creates a quantum state (`State`) in the given `stateSpace` using the provided `value`. This function supports both pure states (`vector`)
and mixed states (`density`), and allows initialization using either amplitudes or Bloch sphere parameters.

# Arguments
- `stateSpace::StateSpace` - The target quantum state space. Must be either `vector` (for pure states) or `density` (for mixed states).
- `value` - Array or structured data describing the initial qubit state(s). Format depends on `stateSpace` and `blochRepresentation`.
- `indexType::IndexType` (optional) - Determines the qubit indexing convention. Defaults to `circuitIndexBigEndian`.
- `blochRepresentation::Bool` (optional) - If `true`, interprets values using Bloch sphere angles ``(θ, ϕ)``. Defaults to `false`.

# Returns
- `State`: A `VectorState` or `DensityState` depending on the chosen `stateSpace`.

# Pure States (`stateSpace` = `vector`)

## Representations

- **Amplitude representation** (default):
  A single-qubit state is defined as: ``|ψ⟩ = α|0⟩ + β|1⟩``
  For multiple qubits, the full state is constructed as: ``|ψ⟩ = |ψ₁⟩ ⊗ |ψ₂⟩ ⊗ ... ⊗ |ψₙ⟩``

- **Bloch representation** (`blochRepresentation = true`):
  Each qubit is defined using Bloch sphere angles: ``|ψ⟩ = cos(θ/2)|0⟩ + sin(θ/2)·exp(i·ϕ)|1⟩``

## Accepted Formats
```julia-repl
createInitialQubitState(vector, [α₁ β₁; ... ; αₙ βₙ])
createInitialQubitState(vector, [αₙ βₙ; ... ; α₁ β₁], indexType=circuitIndexLittleEndian)
createInitialQubitState(vector, [α₁ ... αₙ; β₁ ... βₙ], indexType=byteIndex)

createInitialQubitState(vector, [θ₁ ϕ₁; ... ; θₙ ϕₙ], blochRepresentation=true)
createInitialQubitState(vector, [θₙ ϕₙ; ... ; θ₁ ϕ₁], indexType=circuitIndexLittleEndian, blochRepresentation=true)
createInitialQubitState(vector, [θ₁ ... θₙ; ϕ₁ ... ϕₙ], indexType=byteIndex, blochRepresentation=true)
```

# Mixed States (`stateSpace` = `density`)

## Representations

- **Amplitude representation** (default):
  Mixed states are convex combinations of pure states: ``ρ = p₁|ψ₁⟩⟨ψ₁| + p₂|ψ₂⟩⟨ψ₂| + ... + pₘ|ψₘ⟩⟨ψₘ|``
  Each ``|ψₖ⟩`` is built from individual qubits: ``|ψₖ⟩ = |ψₖ₁⟩ ⊗ |ψₖ₂⟩ ⊗ ... ⊗ |ψₖₙ⟩``
  With each qubit: ``|ψₖⱼ⟩ = αₖⱼ|0⟩ + βₖⱼ|1⟩``

- **Bloch representation** (`blochRepresentation = true`):
  Each qubit in each component state is defined using: ``|ψₖⱼ⟩ = cos(θₖⱼ/2)|0⟩ + sin(θₖⱼ/2)·exp(i·ϕₖⱼ)|1⟩``

## Accepted Formats
```julia-repl
createInitialQubitState(density, [[p₁, [α₁₁ β₁₁; ... ; α₁ₙ β₁ₙ]], ..., [pₘ, [αₘ₁ βₘ₁; ... ; αₘₙ βₘₙ]]])
createInitialQubitState(density, [[p₁, [α₁ₙ β₁ₙ; ... ; α₁₁ β₁₁]], ..., [pₘ, [αₘₙ βₘₙ; ... ; αₘ₁ βₘ₁]]], indexType=circuitIndexLittleEndian)
createInitialQubitState(density, [[p₁, [α₁₁ ... α₁ₙ; β₁₁ ... β₁ₙ]], ..., [pₘ, [αₘ₁ ... αₘₙ; βₘ₁ ... βₘₙ]]], indexType=byteIndex)

createInitialQubitState(density, [[p₁, [θ₁₁ ϕ₁₁; ... ; θ₁ₙ ϕ₁ₙ]], ..., [pₘ, [θₘ₁ ϕₘ₁; ... ; θₘₙ ϕₘₙ]]], blochRepresentation=true)
createInitialQubitState(density, [[p₁, [θ₁ₙ ϕ₁ₙ; ... ; θ₁₁ ϕ₁₁]], ..., [pₘ, [θₘₙ ϕₘₙ; ... ; θₘ₁ ϕₘ₁]]], indexType=circuitIndexLittleEndian, blochRepresentation=true)
createInitialQubitState(density, [[p₁, [θ₁₁ ... θ₁ₙ; ϕ₁₁ ... ϕ₁ₙ]], ..., [pₘ, [θₘ₁ ... θₘₙ; ϕₘ₁ ... ϕₘₙ]]], indexType=byteIndex, blochRepresentation=true)
```

# See Also
- [`IndexType`](@ref)
- [`StateSpace`](@ref)
- [`State`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)

# Examples

**Pure States (vector)**

**Example 1: Big-endian, ``|ψ⟩ = |0⟩⊗|1⟩ → [0; 1; 0; 0]``**
```julia-repl
julia> createInitialQubitState(vector, [1. 0.; 0. 1.])
```
**Example 2: Little-endian, same state**
```julia-repl
julia> createInitialQubitState(vector, [0. 1.; 1. 0.], indexType=circuitIndexLittleEndian)
```
**Example 3: Byte index, ``|ψ⟩ = |0⟩⊗|1⟩⊗|1⟩ → [0; 0; 0; 1; 0; 0; 0; 0]``**
```julia-repl
julia> createInitialQubitState(vector, [1. 0. 0.; 0. 1. 1.], indexType=byteIndex)
```
**Example 4: Byte index for 3 qubits, ``|ψ⟩ = |6⟩ = |1⟩⊗|1⟩⊗|0⟩ → [0; 0; 0; 0; 0; 0; 1; 0]``**
```julia-repl
julia> createInitialQubitState(vector, createByteIndexVector(6, 3), indexType=byteIndex)
```
**Example 5: Big-endian, Bloch representation: ``|ψ⟩ = |+⟩⊗|–⟩⊗|1⟩``**
```julia-repl
julia> createInitialQubitState(vector, [π/2 0.; π/2 π/2], blochRepresentation=true)
```

**Mixed States (density)**

**Example 7: Big-endian, ``ρ = 0.25|+⟩⟨+| + 0.75|–⟩⟨–| → [0.5, -0.25; -0.25, 0.5]``**
```julia-repl
julia> createInitialQubitState(density, [[0.25, [1. 1.]/sqrt(2)], [0.75, [1. -1.]/sqrt(2)]])
```
**Example 8: Same using Bloch representation**
```julia-repl
julia> createInitialQubitState(density, [[0.25, [π/2 0.]], [0.75, [π/2 π]]], blochRepresentation=true)
```

# Notes
The qubit order interpretation depends on indexType. Use this to align with simulation or hardware layout.
All input vectors/matrices must be properly normalized or interpreted as unnormalized (e.g., Bloch angles).
"""
createInitialQubitState(stateSpace::StateSpace, value; indexType = circuitIndexBigEndian, blochRepresentation = false) = dispatchCreateInitialQubitState(Val(stateSpace), value, Val(indexType), Val(blochRepresentation))

dispatchCreateInitialQubitState(::Val{vector}, ab, indexType, blochRepresentation) = VectorState(buildStateVector(ab, indexType, blochRepresentation))
dispatchCreateInitialQubitState(::Val{density}, pab, indexType, blochRepresentation) = DensityState(buildDensityMatrix(pab, indexType, blochRepresentation))

function buildStateVector(ab, indexType, blochRepresentation)
    n = getQubitCount(ab, indexType)
    q = ones(ComplexF64, 2^n, 1)
    for r in 1:(2^n)
        bits = getBitPattern(r, n, indexType)
        q[r] = computeAmplitude(bits, ab, indexType, blochRepresentation)
    end
    return q
end

function buildDensityMatrix(pab, indexType, blochRepresentation)
    n = getQubitCount(pab[1][2], indexType)
    rho = zeros(ComplexF64, 2^n, 2^n)
    for (weight, ab) in pab
        q = buildStateVector(ab, indexType, blochRepresentation)
        rho += weight * q * q'
    end
    return rho
end

function computeAmplitude(bits, ab, indexType, blochRepresentation)
    amp = one(ComplexF64)
    for c in eachindex(bits)
        a, b = getAmplitudePair(ab, c, indexType)
        amp *= bits[c] == 0 ? firstComponent(a, b, blochRepresentation) : secondComponent(a, b, blochRepresentation)
    end
    return amp
end

firstComponent(a, _, ::Val{false}) = a
firstComponent(a, _, ::Val{true}) = cos(a / 2)
secondComponent(_, b, ::Val{false}) = b
secondComponent(a, b, ::Val{true}) = exp(im * b) * sin(a / 2)

getAmplitudePair(ab, c, ::Val{byteIndex}) = (ab[1, c], ab[2, c])
getBitPattern(r, n, ::Val{byteIndex}) = reverse(digits(r - 1, base = 2, pad = n)')
getQubitCount(ab, ::Val{byteIndex}) = size(ab, 2)

getAmplitudePair(ab, c, ::Val{circuitIndexBigEndian}) = (ab[c, 1], ab[c, 2])
getBitPattern(r, n, ::Val{circuitIndexBigEndian}) = reverse(digits(r - 1, base = 2, pad = n)')
getQubitCount(ab, ::Val{circuitIndexBigEndian}) = size(ab, 1)

getAmplitudePair(ab, c, ::Val{circuitIndexLittleEndian}) = (ab[c, 1], ab[c, 2])
getBitPattern(r, n, ::Val{circuitIndexLittleEndian}) = digits(r - 1, base = 2, pad = n)'
getQubitCount(ab, ::Val{circuitIndexLittleEndian}) = size(ab, 1)

"""
    createDoubleQubitWernerDensityState(
		p::Real; 
		indexType::IndexType = circuitIndexBigEndian
	) -> DensityState

Creates a **Werner density state** for two qubits with mixing parameter ``p``.

The Werner density state is defined as a mixture of a maximally entangled singlet state and the maximally mixed state: ``ρ = p·|Ψ⁻⟩⟨Ψ⁻| + (1 - p)·I/4``

where ``|Ψ⁻⟩ = (|01⟩ - |10⟩) / √2`` is the **singlet Bell state**, and ``I`` is the 4×4 identity matrix.

# Arguments
- `p::Real` - Mixing probability ``p``. Must be in the range ``[0, 1]``.
  - ``p = 1``: Fully entangled (pure singlet state).
  - ``p = 0``: Completely mixed state (maximally impure).
- `indexType::IndexType` (optional) - Determines qubit indexing convention. Defaults to `circuitIndexBigEndian`.

# Returns
- `DensityState` - A two-qubit mixed quantum state ``ρ`` represented as a `DensityState` object.

# Properties
- **Bell inequality** is violated if ``p > 1/√2``.
- **Positive partial transpose (PPT)** criterion indicates entanglement if ``p > 1/3``.

# See Also
- [`IndexType`](@ref)
- [`DensityState`](@ref)

# Examples

**Example 1**: Completely mixed state ``p = 0``
```julia-repl
julia> createDoubleQubitWernerDensityState(0).rho
4×4 Matrix{ComplexF64}:
 0.25+0.0im   0.0+0.0im   0.0+0.0im   0.0+0.0im
  0.0+0.0im  0.25+0.0im   0.0+0.0im   0.0+0.0im
  0.0+0.0im   0.0+0.0im  0.25+0.0im   0.0+0.0im
  0.0+0.0im   0.0+0.0im   0.0+0.0im  0.25+0.0im
```

**Example 2**: Pure singlet state ``p = 1``
```julia-repl
julia> createDoubleQubitWernerDensityState(1).rho
4×4 Matrix{ComplexF64}:
 0.0+0.0im   0.0+0.0im   0.0+0.0im  0.0+0.0im
 0.0+0.0im   0.5+0.0im  -0.5+0.0im  0.0+0.0im
 0.0+0.0im  -0.5+0.0im   0.5+0.0im  0.0+0.0im
 0.0+0.0im   0.0+0.0im   0.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitWernerDensityState(p::Real; indexType::IndexType = circuitIndexBigEndian)
	# rho=p*|singlet><singlet|+0.25*(1-p)*identity, here we choose |singlet>=(|01>-|10>)/sqrt(2)
	# 1=pure entangled, 0=impure, Bell's inequality is satisfied if p<1/sqrt(2), 
	# PPT condition is satisfied if p<1/3
	return DensityState(p * [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0] / 2 + (1-p) * [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1] / 4)
end

#################################################################################################################################

mutable struct Cache
	cV1::ComplexMatrix
	cV2::ComplexMatrix
	cM1::ComplexMatrix
	cM2::ComplexMatrix
	outcome::Int
end

function applyOperationOnQubitState(unitaryOperation::UnitaryOperation, vectorState::VectorState, cache::Cache)
	return VectorState(unitaryOperation.U * vectorState.q)
end

function applyOperationOnQubitState(unitaryOperation::UnitaryOperation, densityState::DensityState, cache::Cache)
    mul!(cache.cM1, densityState.rho, (unitaryOperation.U'))
	return DensityState(unitaryOperation.U * cache.cM1)
end

function applyOperationOnQubitState(measureOperation::MeasureOperation, vectorState::VectorState, cache::Cache)
	prob = zeros(length(measureOperation.krausOperators))
	mul!(cache.cM1, vectorState.q, vectorState.q')
	for m in axes(measureOperation.krausOperators, 1)
		mul!(cache.cM2, measureOperation.measurementOperators[m], cache.cM1)
		prob[m] = real(tr(cache.cM2))
	end
	cache.outcome = sampleMeasurementOutcome(prob)
	mul!(cache.cV1, vectorState.q, 1 / sqrt(prob[cache.outcome+1]))
	return VectorState(measureOperation.measurementOperators[cache.outcome+1] * cache.cV1; measured = Measured(prob, cache.outcome, measureOperation.labels))
end

function applyOperationOnQubitState(measureOperation::MeasureOperation, densityState::DensityState, cache::Cache)
	prob = zeros(length(measureOperation.krausOperators))
	for m in axes(measureOperation.krausOperators, 1)
		mul!(cache.cM2, measureOperation.measurementOperators[m], densityState.rho)
		prob[m] = real(tr(cache.cM2))
	end
	cache.outcome = sampleMeasurementOutcome(prob)
	mul!(cache.cM1, measureOperation.krausOperators[cache.outcome+1].E', 1 / prob[cache.outcome+1])
	mul!(cache.cM2, densityState.rho, cache.cM1)
	return DensityState(measureOperation.krausOperators[cache.outcome+1].E * cache.cM2; measured = Measured(prob, cache.outcome, measureOperation.labels))
end

function applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, vectorState::VectorState, cache::Cache)
	throw("Combination of a MeasureAndForgetOperation with a VectorState is not possible...")
end

function applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, densityState::DensityState, cache::Cache)
	prob = zeros(length(measureAndForgetOperation.krausOperators))
	rhoa = zeros(size(densityState.rho))
	for m in axes(measureAndForgetOperation.krausOperators, 1)
		mul!(cache.cM1, measureAndForgetOperation.measurementOperators[m], densityState.rho)
		prob[m] = real(tr(cache.cM1))
		mul!(cache.cM1, densityState.rho, measureAndForgetOperation.krausOperators[m].E')
		mul!(cache.cM2, measureAndForgetOperation.krausOperators[m].E, cache.cM1)
		rhoa .+= cache.cM2
	end
	return DensityState(rhoa; measured = Measured(prob, nothing, measureAndForgetOperation.labels))
end

function sampleMeasurementOutcome(probabilities::Vector{Float64})
	p = rand()
	cumulative = cumsum(probabilities, dims=1)
	indices = findmin(findall(x -> p <= x, cumulative))
	return first(indices)-1
end

function applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, vectorState::VectorState, cache::Cache)
	throw("Combination of a QuantumChannelOperation with a VectorState is not possible...")
end

function applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, densityState::DensityState, cache::Cache)
	rhoa = zeros(size(densityState.rho))
	for krausOperator in quantumChannelOperation.krausOperators
		mul!(cache.cM1, densityState.rho, krausOperator.E')
		mul!(cache.cM2, krausOperator.E, cache.cM1)
		rhoa .+= cache.cM2
	end
	return DensityState(rhoa)
end

""" 
    runQuantumProgram(
		quantumProgram::QuantumProgram, 
		initialQubitState::S, 
		numberOfShots::Int
	) -> QuantumOutput{Z, I, S} where {Z, I, S <: State}

Executes a quantum program multiple times using a specified initial quantum state.

This function runs the provided `quantumProgram` starting from the `initialQubitState` and simulates its execution `numberOfShots` times. The result includes the full quantum state evolution across all operations for each shot.

This function is **compile-time specialized** on type parameter `S <: State` (that forces the struct to hold type-stable quantum state objects), ensuring optimized and branch-free code generation.

# Type parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.
- `S <: State`: The concrete type of the quantum states being stored, e.g., `VectorState` or `DensityState`.

# Arguments
- `quantumProgram::QuantumProgram` - A compiled quantum program containing a sequence of quantum operations.
- `initialQubitState::S` - The starting state for all shots. This can be a `VectorState` (pure state) or `DensityState` (mixed state).
- `numberOfShots::Int` - Number of repetitions (shots) to run the program. Each shot evolves independently.

# Supported State Types
- `VectorState` - For pure states.
- `DensityState` - For mixed states.

# Supported Operations
- `UnitaryOperation` - Applies a reversible gate.
- `MeasureOperation` - Simulates measurement with revealing outcome by collapsing the state to one of the possible outcomes.
- `MeasureAndForgetOperation` - Simulates measurement without revealing outcome by creating a mixed state over all possible outcomes.
- `QuantumChannelOperation` - Applies a general (possibly non-unitary) quantum channel.

# Returns
- `QuantumOutput{Z, I, S}` - Contains all intermediate states across all operations and shots. The `.output` field is a `Matrix{S}` of size `(numberOfStates, numberOfShots)`.

# Output Layout
If the program has `n` operations, the output is an `(n+1) × numberOfShots` matrix, where:
- Row `1` contains the initial state.
- Row `2` through `n+1` contain the resulting states after each operation.
- Each column corresponds to a different shot.

# Notes
- Each shot is processed independently; the program is re-applied from scratch for each shot.
- The function applies the sequence of operations to the initial state and stores the intermediate results for analysis or visualization.
- Measurement outcomes and post-measurement states are captured in the final states.

# Quantum state evolution 

- Under a unitary operation
    - Let U be the unitary matrix defining the unitary operation.
        - **Quantum vector state**:  
            - A pure state ``|ψ⟩`` evolves as: ``|ψ⟩ → U·|ψ⟩``
        - **Quantum density state**:
            - A density operator ``ρ`` evolves as: ``ρ → U·ρ·U^†``

- Under a generalized measurement (POVM)
    - Let ``{Mₖ}`` be measurement operators satisfying ``∑ₖ Mₖ^†·Mₖ = I``. Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
        - **Quantum vector state**:
            - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·|ψ⟩⟨ψ|)``
            - Post-measurement state for outcome ``k``: ``|ψ⟩ → Mₖ·|ψ⟩/√pₖ``
        - **Quantum density state**:
            - Probability of outcome ``k``: ``pₖ = Tr(Mₖ·ρ)``
            - Post-measurement state (outcome ``k`` observed): ``ρ → Eₖ·ρ·Eₖ^†/pₖ``
            - Post-measurement state (outcome not observed): ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

- Under a projective measurement (PVM)
    - A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``. Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`. Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

- Under a quantum channel
    - Let ``{Eₖ}`` be the Kraus operators of the quantum channel.
        - **Quantum density state**:
            - A density operator ``ρ`` evolves as: ``ρ → ∑ₖ Eₖ·ρ·Eₖ^†``

# See Also
- [`QuantumProgram`](@ref)
- [`QuantumOutput`](@ref)
- [`State`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)
- [`UnitaryOperation`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)

# Examples

## Example 1: Single-qubit vector initial state ``|0⟩``, Hadamard gate, then measurement. Take 5 shots and extract the last state after the measurement.
```julia-repl
julia> qc=createQuantumCircuit(1)
julia> hGate!(qc,1)
julia> measureGate!(qc,[1])
julia> qp=compileQuantumCircuit(qc)
julia> iqs=createInitialQubitState(vector, [1. 0.;])
julia> qo=runQuantumProgram(qp, iqs, 5)
julia> getState(qo, 3)
5-element Vector{VectorState}:
 VectorState(ComplexF64[1.0 + 0.0im; 0.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
 VectorState(ComplexF64[0.0 + 0.0im; 1.0 + 0.0im;;], Measured([0.5000000000000001, 0.4999999999999999], 1, ["0", "1"]))
```

## Example 2: Same circuit with a zero mixed (density) input state. Take 5 shots and extract the last state after the measurement.
```julia-repl
julia> qc=createQuantumCircuit(1)
julia> hGate!(qc,1)
julia> measureGate!(qc,[1])
julia> qp=compileQuantumCircuit(qc)
julia> iqs=createInitialQubitState(density, [[1.00,[1. 0.;]]])
julia> qo=runQuantumProgram(qp, iqs, 5)
julia> getState(qo, 3)
5-element Vector{DensityState}:
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
 DensityState(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], 0.0, 1.0, 1.0, Measured([0.5000000000000001, 0.4999999999999999], 0, ["0", "1"]))
```
"""
function runQuantumProgram(quantumProgram::QuantumProgram, initialQubitState::S, numberOfShots::Int) where {S <: State}
	N = 2^quantumProgram.numberOfQubits
	cache = Cache(ComplexMatrix(undef, N, 1), ComplexMatrix(undef, N, 1), ComplexMatrix(undef, N, N), ComplexMatrix(undef, N, N), 0)
	numberOfStates = length(quantumProgram.program)+1

    # Create a type-stable matrix with concrete type S. We use `undef` because we will fill it in the loop.
    output = Matrix{S}(undef, numberOfStates, numberOfShots)

    # Initialize the first state for all shots (using broadcasting (.=) for VectorState and DensityState).
    output[1, :] .= initialQubitState

    # Loop in cache-friendly order (over columns inside loopOverShots!)
	for (m, operation) in enumerate(quantumProgram.program)
		# call _loopOverShots with type-stable operation::O where {O <: Operation}
		loopOverShots!(operation, output, m, numberOfShots, cache)
    end

	return QuantumOutput(quantumProgram.settings, quantumProgram.numberOfQubits, numberOfStates, numberOfShots, output)
end

function loopOverShots!(operation::O, output::Matrix{S}, m::Int, numberOfShots::Int, cache::Cache) where {O <: Operation, S <: State}
    @inbounds @simd for k in 1:numberOfShots
        output[m+1, k] = applyOperationOnQubitState(operation, output[m, k], cache)
    end
end

""" 
    calculateEntropyAndPurity(
		rho; 
		accuracyCheckForPurity = 10 * eps(1.0)
	) -> (entropy, purity, purityNorm, purityString)

Computes the **von Neumann entropy** and **purity** of a quantum density matrix ``ρ``.

This function analyzes the degree of *mixedness* of a quantum state. The entropy quantifies the uncertainty (disorder), and the purity characterizes how close the state is to being a pure state. Results include a normalized purity and a human-readable description of the state’s purity level.

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `accuracyCheckForPurity::Float64` (optional) - Numerical tolerance to determine whether the purity is effectively 1 or `1/n`. Defaults to `10 * eps(1.0)`.

# Returns
- `entropy::Float64` - Von Neumann entropy, computed as ``−Tr(ρ·log(ρ))``. Zero for pure states.
- `purity::Float64` - Raw purity value: ``Tr(ρ²)``, where:
  - ``1/n`` indicates a maximally mixed state,
  - ``1`` indicates a pure state.
- `purityNorm::Float64` - Normalized purity in range ``[0, 1]``, where:
  - ``0`` indicates a maximally mixed state,
  - ``1`` indicates a pure state.
- `purityString::String` - Descriptive label of purity:
  - `"pure"` if ``ρ`` is approximately pure,
  - `"impure, totally mixed"` if ``ρ`` is close to ``(1/n)·I``,
  - `"impure, partially mixed"` otherwise.

# Purity Range
- Pure state: purity = ``1``
- Maximally mixed state: purity = ``1 / n``
- Normalized purity = ``(purity·n - 1) / (n - 1)``

# Notes
- Entropy uses eigenvalue decomposition: entropy = ``−∑ᵢ λᵢ·log(λᵢ)``, skipping ``λᵢ ≈ 0`` or ``1``.
- The function assumes ``ρ`` is a valid density matrix (Hermitian, positive semi-definite, unit trace).
- The purityString classification is tolerant to floating-point precision via accuracyCheckForPurity.

# See Also
- [`DensityState`](@ref)

# Examples

**Example 1**: Totally mixed density state
```julia-repl
julia> (entropy, purity, purityNorm, purityString)=calculateEntropyAndPurity([0.5 0.0; 0.0 0.5])
(0.6931471805599453, 0.5, 0.0, "impure, totally mixed")
```

**Example 2**: Pure density state
```julia-repl
julia> (entropy, purity, purityNorm, purityString)=calculateEntropyAndPurity([0.5 0.5; 0.5 0.5])
(0, 1.0, 1.0, "pure")
```
"""
function calculateEntropyAndPurity(rho; accuracyCheckForPurity = 10*eps(1.0))
	n = size(rho, 1)
	# accuracyCheckForPurity = 1e-6
	entropy = 0
	purity = 0
	L = eigvals(rho)
	for i in 1:n
		purity = purity + L[i]^2 # purity=Tr(rho^2)
		if ((L[i] != 0) && (L[i] != 1))
			entropy = entropy - L[i] * log(abs(L[i])) # entropy=-Tr(rho*ln(rho))
			# use limit x*log(x) --> 0 (x --> 0)
			# use limit x*log(x) --> 0 (x --> 1)
		end
	end
	# for finite system, entropy quantifies departure of system from a pure state. 
	# it codifies degree of mixing of given finite state. for pure state entropy==0
	# purity is in the range [1/n,...,1] where [1/n=impure totally mixed, impure partially mixed, 1=pure]
	purityNorm = (purity * n - 1) / (n - 1) # purityNorm(purity)=0+((1-0)/(1-1/n))*(purity-1/n)
	if abs(purity - (1 / n)) < accuracyCheckForPurity
		purityString = "impure, totally mixed"
	elseif abs(purity - 1) < accuracyCheckForPurity
		purityString = "pure"
	else
		purityString = "impure, partially mixed"
	end
	return real(entropy), real(purity), real(purityNorm), purityString
end

# """ 
#     calculateExpectationValueOfQubitState(q) -> Vector{Float64}

# Computes the expectation value (probability distribution) of a qubit state ``|q⟩``.

# This function returns a vector where each element represents the probability of measuring the qubit in the corresponding basis state. The probabilities are computed as the squared magnitudes of the complex amplitudes in ``|q⟩``.

# # Arguments
# - `q::AbstractVector` - A normalized qubit state vector ``|q⟩`` (typically of complex numbers).

# # Returns
# - `Vector{Float64}` - A vector of probabilities corresponding to each basis state.

# # Example
# ```julia-repl
# julia> calculateExpectationValueOfQubitState([1; 1]/sqrt(2))
# 2-element Vector{Float64}:
#  0.4999999999999999
#  0.4999999999999999
# ```
# """
function calculateExpectationValueOfQubitState(q)
	return abs.(q).^2
end

"""
    calculateCorrelation(
		x1::AbstractVector, 
		x2::AbstractVector
	) -> Float64

Compute the Pearson correlation coefficient between two vectors `x1` and `x2`.

This function measures the linear correlation between `x1` and `x2`. A value of ``1.0`` indicates perfect positive linear correlation, ``-1.0`` indicates perfect negative linear correlation, and ``0.0`` indicates no linear correlation.

# Arguments
- `x1::AbstractVector` - The first data vector.
- `x2::AbstractVector` - The second data vector. Must be the same length as `x1`.

# Returns
- `Float64` - The Pearson correlation coefficient between `x1` and `x2`.

# Examples
```julia-repl
julia> calculateCorrelation([1 2 3], [3 2 1])
-1.0
```
"""
function calculateCorrelation(x1, x2)
	xa1 = x1.-sum(x1)/length(x1)
	xa2 = x2.-sum(x2)/length(x2)
	sdxa12 = sqrt(sum(xa1.^2)*sum(xa2.^2))
	correlation = sum(xa1.*xa2)/sdxa12
	return correlation
end

"""
    @enum UseFromMeasurement useProbability useOutcome

Specify the source of data for computing expectation values from a quantum measurement.

# Variants
- `useProbability` - Use the theoretical probability distribution from the quantum measurement.
- `useOutcome` - Use the actual sampled outcomes (e.g., from repeated measurements or shots).

This enum is used to control whether to compute expectation values based on ideal probabilities or empirical results.
"""
@enum UseFromMeasurement useProbability useOutcome

"""
    calculateExpectationValueOfProductOfMeasuredQubits(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		useFromMeasurement::UseFromMeasurement
	) -> Float64

Computes the expectation value of the **product** of measured qubits from the given quantum simulation output.

The expectation value is computed by checking the parity of the measurement outcomes:
- If the number of `1`s is even, contribute `+1`.
- If the number of `1`s is odd, contribute `-1`.

# Parity Evaluation Example (2-qubit system)

| Bitstring | Number of 1s | Parity | Contribution |
| --------- | ------------ | ------ | ------------ |
| `00`      | 0            | Even   | `+1`         |
| `01`      | 1            | Odd    | `-1`         |
| `10`      | 1            | Odd    | `-1`         |
| `11`      | 2            | Even   | `+1`         |

# Arguments
- `quantumOutput::QuantumOutput` - The result of running a quantum program.
- `stateId::Int` - The index of the quantum state to extract results from.
- `useFromMeasurement::UseFromMeasurement` - Determines whether to use probabilities (`useProbability`) or actual measurement outcomes (`useOutcome`) for the computation.

# Returns
- `Float64` - The expectation value of the product of the measured qubit values.

# See Also
- [`QuantumOutput`](@ref)
- [`UseFromMeasurement`](@ref)

# Examples

**Uncorrelated 2-qubit superposition state**

Each qubit is put in superposition using a Hadamard gate: ``|ψ⟩ = (|0⟩ + |1⟩)/√2 ⊗ (|0⟩ + |1⟩)/√2 = (1/2)(|00⟩ + |01⟩ + |10⟩ + |11⟩)``.

All bitstrings occur with equal probability ``0.25``. Expectation value is:

| Bitstring | Number of 1s | Parity | Contribution | Probability | Weighted Contribution |
| --------- | ------------ | ------ | ------------ | ----------- | --------------------- |
| `00`      | 0            | Even   | `+1`         | 0.25        | `+0.25`               |
| `01`      | 1            | Odd    | `-1`         | 0.25        | `-0.25`               |
| `10`      | 1            | Odd    | `-1`         | 0.25        | `-0.25`               |
| `11`      | 2            | Even   | `+1`         | 0.25        | `+0.25`               |

Final Expectation Value: ``0.25 − 0.25 − 0.25 + 0.25 = 0.0``

```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> hGate!(qc, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> corr = calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
0.0
```

**Correlated 2-qubit Bell state**

We prepare an entangled Bell state: ``|ψ⟩ = (|00⟩ + |11⟩) / √2``

Only `00` and `11` occur, each with probability ``0.5``. Expectation value is:

| Bitstring | Number of 1s | Parity | Contribution | Probability | Weighted Contribution |
| --------- | ------------ | ------ | ------------ | ----------- | --------------------- |
| `00`      | 0            | Even   | `+1`         | 0.5         | `+0.5`                |
| `01`      | 1            | Odd    | `-1`         | 0.0         | ` 0.0`                |
| `10`      | 1            | Odd    | `-1`         | 0.0         | ` 0.0`                |
| `11`      | 2            | Even   | `+1`         | 0.5         | `+0.5`                |

Final Expectation Value: ``0.5 + 0.0 + 0.0 + 0.5 = 1.0``

```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps = true)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> corr = calculateExpectationValueOfProductOfMeasuredQubits(qo, 3, useProbability)
1.0
```
"""
function calculateExpectationValueOfProductOfMeasuredQubits(quantumOutput::QuantumOutput, stateId::Int, useFromMeasurement::UseFromMeasurement)
    nrMeasuredBits = Int(log2(length(quantumOutput.output[stateId].measured.probability)))
    expVal = 0.
    if useFromMeasurement == useProbability
        for m in 1:2^nrMeasuredBits
            bin = reverse(digits(m-1, base = 2, pad = nrMeasuredBits)')
            if iseven(sum(bin))
                expVal += quantumOutput.output[stateId].measured.probability[m]
            else
                expVal += -quantumOutput.output[stateId].measured.probability[m]
            end
        end
    else # useOutcome
        nrShots = size(quantumOutput.output, 2)
        for k in 1:nrShots
            bin = reverse(digits(quantumOutput.output[stateId, k].measured.outcome, base = 2, pad = nrMeasuredBits)')
            if iseven(sum(bin))
                expVal += 1.
            else
                expVal += -1.
            end
        end
        expVal = expVal/nrShots
    end
    return expVal
end

""" 
    generateKrausOperatorsForPVMMeasurement(
        numberOfQubits::Int,
        qubits::Qubits,
        sigmas::Sigmas;
        zeroBasedNumbering::Bool = false,
        indexType::IndexType = circuitIndexBigEndian
    ) -> KrausOperators

Generates Kraus operators for a **projective measurement (PVM)** on a subset of qubits in a quantum circuit.

This function constructs the Kraus operators corresponding to a measurement performed on the specified `qubits` within a system of `numberOfQubits` qubits.

The measurement directions are defined by the list of operators `sigmas`. If `sigmas` is empty, the default ``σ_z`` is used for each qubit.

# Arguments
- `numberOfQubits::Int` - Total number of qubits in the system.
- `qubits::Qubits` - Indices of the qubits to be measured, interpreted using the settings for `zeroBasedNumbering` and `indexType`.
- `sigmas::Sigmas` - Associated measurement direction operators for each qubit. If `sigmas` is empty (`[]`), ``σ_z`` is used for all specified qubits.
- `zeroBasedNumbering::Bool` (optional): 
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1. Defaults to `false`.
- `indexType::IndexType` (optional) - Determines qubit indexing convention. Defaults to `circuitIndexBigEndian`.

# Returns
- `KrausOperators` — List of Kraus operators corresponding to the **projective measurement (PVM)**.

# How **PVM** is related to **POVM**
- A **POVM** is a generalized measurement defined by measurement operators ``{Mₖ}`` satisfying ``∑ₖ Mₖ^†·Mₖ = I``. Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
- A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``. Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`. Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

# See Also
- [`IndexType`](@ref)
- [`Qubits`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
We prepare Kraus operators for a **projective measurement (PVM)** on a single qubit, using the measurement operator ``σ_z``, which corresponds to a measurement in the computational basis ``{|0⟩, |1⟩}``.
This results in the following Kraus operators: ``E_0 = |0⟩⟨0| = [1, 0; 0, 0]``, ``E_1 = |1⟩⟨1| = [0, 0; 0, 1]``
```julia-repl
julia> generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()])
2-element Vector{Any}:
 ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
```
"""
function generateKrausOperatorsForPVMMeasurement(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas; zeroBasedNumbering = false, indexType::IndexType = circuitIndexBigEndian)
	# assume bits conforms to Settings
	settings = Settings(zeroBasedNumbering, indexType)
	qubits = convertToByteIndex(settings, numberOfQubits, qubits)
	sigmas = sigmas == [] ? [sigmaZ() for _ in qubits] : sigmas
	krausOperators = constructKrausOperatorsForPVMMeasurement(numberOfQubits, qubits, sigmas)
	return krausOperators
end

""" 
    generateKrausOperatorsForDepolarizingChannel(
		p::Float64
	) -> KrausOperators

Constructs the Kraus operators for a single-qubit **depolarizing quantum channel**.  
This channel applies a random Pauli error ``σ_x``, ``σ_y`` or ``σ_z`` with total probability ``p``, and leaves the state unchanged with probability ``1 - p``.
The depolarizing channel is defined as:

```math
ρ → (1-p)·ρ + (p/3)·(σ_x·ρ·σ_x+σ_y·ρ·σ_y+σ_z·ρ·σ_z)
```

where ``ρ`` is the input density matrix, ``σ_x``, ``σ_y`` and ``σ_z`` are the Pauli matrices and ``p`` is the error probability.
The corresponding Kraus operators are: ``M_I = √(1-p)·I``, ``M_X = √(p/3)·σ_x``, ``M_Y = √(p/3)·σ_y`` and ``M_Z = √(p/3)·σ_z``.

# Effects
- Both populations driven toward ``½``
- Off-diagonals shrink by ``1 - 4p/3``

# Arguments
- `p::Float64` - The probability that a qubit undergoes an error.

# Returns
- `KrausOperators` — List of Kraus operators ``M_I, M_X, M_Y, M_Z`` corresponding to the **depolarizing quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying a depolarizing channel with error probability `p = 0.75` on a single qubit in superposition:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForDepolarizingChannel(0.75))
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> # First state: pure superposition (|0⟩ + |1⟩)/√2
julia> round.(real.(qo.output[1].rho), digits=3)
[0.5 0.5; 0.5 0.5]
julia> qo.output[1].entropy
-9.056320547599769e-31
julia> qo.output[1].purityValue
1.0
julia> qo.output[1].purityNorm
1.0

julia> # Second state: completely mixed state (Id/2) due to p = 0.75
julia> round.(real.(qo.output[2].rho), digits=3)
[0.5 0.0; 0.0 0.5]
julia> qo.output[2].entropy
0.6931471805599453
julia> qo.output[2].purityValue
0.5
julia> qo.output[2].purityNorm
0.0
```
"""
function generateKrausOperatorsForDepolarizingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(sqrt(1-p)*ComplexMatrix([1. 0; 0. 1.]), "MI"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliX(), "MX"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliY(), "MY"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliZ(), "MZ"))
	return krausOperators
end

""" 
    generateKrausOperatorsForPhaseDampingChannel(
        p::Float64
    ) -> KrausOperators

Constructs the Kraus operators for a single-qubit phase damping quantum channel.
This channel models the loss of quantum coherence (dephasing) without energy dissipation.
With probability ``p``, the off-diagonal elements of the density matrix decay toward zero, while populations remain unchanged.
The phase damping channel acts on a density matrix ``ρ`` as:

```math
ρ = [ρ_{00}, ρ_{01}; ρ_{10}, ρ_{11}] → [ρ_{00}, (1−p)⋅ρ_{01}; (1−p)⋅ρ_{10}, ρ_{11}]
```

where
- the diagonal elements (populations) are unaffected,,
- the off-diagonal elements (coherences) decay linearly with factor ``1-p``.

The corresponding Kraus operators are: ``M_0 = √(1-p)·I`` (induces no decoherence), ``M_1 = √p·|0⟩⟨0|`` (induces decoherence) and ``M_2 = √p·|1⟩⟨1|`` (induces decoherence).

# Effects
- Populations unchanged
- Off-diagonals shrink linearly by ``1-p``

# Arguments
- `p::Float64` - The probability that coherence between ``|0⟩`` and ``|1⟩`` is lost.

# Returns
- `KrausOperators` — List of Kraus operators ``M_0, M_1, M_2`` corresponding to the **phase damping quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying a phase damping channel with p = 0.5 on a qubit in superposition:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPhaseDampingChannel(0.5))
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> round.(real.(qo.output[2].rho), digits=3)
[0.5 0.25; 0.25 0.5]
```
"""
function generateKrausOperatorsForPhaseDampingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(sqrt(1-p)*ComplexMatrix([1. 0; 0. 1.]), "M0"))
	push!(krausOperators, KrausOperator(sqrt(p)*ComplexMatrix([1. 0; 0. 0.]), "M1"))
	push!(krausOperators, KrausOperator(sqrt(p)*ComplexMatrix([0. 0; 0. 1.]), "M2"))
	return krausOperators
end

""" 
    generateKrausOperatorsForAmplitudeDampingChannel(
        p::Float64
    ) -> KrausOperators

Constructs the Kraus operators for a single-qubit amplitude damping quantum channel.
This channel models energy dissipation, where the excited state ``|1⟩`` decays to the ground state ``|0⟩`` with probability ``p``.
The amplitude damping channel acts on a density matrix ``ρ`` as:

```math
ρ = [ρ_{00}, ρ_{01}; ρ_{10}, ρ_{11}] → [ρ_{00}+p⋅ρ_{11}, √(1−p)⋅ρ_{01}; √(1−p)⋅ρ_{10}, (1−p)⋅ρ_{11}]
```

where
- the population of ``|1⟩`` decays into ``|0⟩`` with probability ``p``,
- off-diagonal coherence terms shrink by ``√(1-p)``,
- population of ``|0⟩`` increases accordingly, preserving trace.

The corresponding Kraus operators are: ``M_0 = |0⟩⟨0| + √(1-p)⋅|1⟩⟨1|`` (induces no decay) and ``M_1 = √p⋅|0⟩⟨1|`` (induces decay ``|1⟩ → |0⟩``).

# Effects
- ``ρ_{11} → (1-p)ρ_{11}``, population leaks to ``ρ_{00}``
- Off-diagonals shrink by ``√(1-p)``

# Arguments
- `p::Float64` - The probability of spontaneous decay ``|1⟩ → |0⟩``.

# Returns
- `KrausOperators` — List of Kraus operators ``M_0, M_1`` corresponding to the **amplitude damping quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying an amplitude damping channel with p = 0.3 on a qubit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForAmplitudeDampingChannel(0.3))
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> round.(real.(qo.output[2].rho), digits=3)
[0.65 0.418; 0.418 0.35]
```
"""
function generateKrausOperatorsForAmplitudeDampingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(ComplexMatrix([1. 0; 0. sqrt(1-p)]), "M0"))
	push!(krausOperators, KrausOperator(ComplexMatrix([0. sqrt(p); 0. 0.]), "M1"))
	return krausOperators
end

""" 
    partialTrace(
		rho::Matrix{<:Complex}; 
		traceIndex::Int, 
		zeroBasedNumbering = false, 
		indexType = circuitIndexBigEndian
	) -> Matrix{<:Complex}

Computes the partial trace of a given density matrix ``ρ`` over a specified subsystem index `traceIndex`. 
This operation removes (traces out) one qubit from a larger system, resulting in a reduced density matrix for the remaining subsystem(s).

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `traceIndex::Int` — The qubit index to trace out, interpreted using the settings for `zeroBasedNumbering` and `indexType`.
- `zeroBasedNumbering::Bool` (optional) — Whether `traceIndex` uses zero-based indexing. Defaults to `false`.
- `indexType::IndexType` (optional) — Interpretation of `traceIndex` based on circuit/byte ordering. Defaults to `circuitIndexBigEndian`.

# Returns
- `Matrix{<:Complex}` — Reduced density matrix.

# See Also
- [`IndexType`](@ref)

# Example
Trace out the first qubit from the Bell state ``|ψ⟩ = (|00⟩ + |11⟩) / √2``:
```julia-repl
julia> partialTrace([1 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 1]/2., traceIndex = 1)
2×2 Matrix{Float64}:
 0.5  0.0
 0.0  0.5
```
"""
function partialTrace(rho; traceIndex::Int, zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
	numberOfQubits = Int(log2(size(rho, 1)))

	# assume traceIndex conforms to Settings
	traceIndex = convertToByteIndex(Settings(zeroBasedNumbering, indexType), numberOfQubits, traceIndex)

	list0 = zeros(Int, 2^(numberOfQubits-1), 1)
	list1 = zeros(Int, 2^(numberOfQubits-1), 1)
	p = 1
	q = 1
	for m in 1:2^numberOfQubits
		bin = reverse(digits(m-1, base = 2, pad = numberOfQubits)')
		if bin[traceIndex] == 0
			list0[p] = m
			p += 1
		else
			list1[q] = m
			q += 1
		end
	end

	T0 = zeros(Int, 2^(numberOfQubits-1), 2^numberOfQubits)
	T1 = zeros(Int, 2^(numberOfQubits-1), 2^numberOfQubits)
	for m in 1:2^(numberOfQubits-1)
		T0[m, list0[m]] = 1
		T1[m, list1[m]] = 1
	end

	return T0 * rho * T0' + T1 * rho * T1'
end

""" 
    fidelity(
		rho::Matrix{<:Complex}, 
		sigma::Matrix{<:Complex}
	) -> Float64

Compute the **quantum fidelity** between two density matrices ``ρ`` and ``σ``.

Fidelity is a measure of similarity between quantum states, defined as: ``F(ρ, σ) = Tr[√ρ·σ·√ρ]``

# Arguments
- `rho::Matrix{<:Complex}` - Density matrix ``ρ`` representing a quantum state (Hermitian, positive semi-definite, trace 1).
- `sigma::Matrix{<:Complex}` - Another density matrix ``σ`` representing a quantum state of the same size as ``ρ``.

# Returns
- `Float64` - The fidelity between ``ρ`` and ``σ``, a value between 0 and 1

# See Also
- [`DensityState`](@ref)

# Example
Compute the fidelity between the pure states ``ρ = |00⟩⟨00|`` and ``σ = |11⟩⟨11|``:
```julia-repl
julia> fidelity([1 0; 0 0], [0 0; 0 1])
0.0
```
"""
function fidelity(rho, sigma)
	return tr(sqrt(rho) * sigma * sqrt(rho))
end

###########################

""" 
    probeMeasureOutcome(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		aTitle::String
	) -> Nothing

Visualizes the **measured outcomes** of a quantum circuit run as a histogram.

This function extracts the measurement results from a specific qubit register or composite system (indexed by `stateId`) within `quantumOutput`, then plots the frequency of each outcome over all measurement shots.
The histogram is labeled with the provided title `aTitle`.

# Arguments
- `quantumOutput::QuantumOutput` - The result object returned by `runQuantumProgram`, containing the outcomes of quantum measurements.
- `stateId::Int` - The ID of the quantum register or subsystem whose measurement results are to be plotted. This will be internally converted to 1-based indexing if necessary.
- `aTitle::String` - The title to display on the histogram plot.

# Returns
- `Nothing` - Makes the plot.

# Notes
- This function assumes the quantum output is structured with measurement results accessible via `quantumOutput.output[stateId, shot]`.
- The x-axis labels are extracted from the measurement label metadata (`measured.labels`) of the first shot.
- The y-axis shows the normalized frequency of each measured bitstring combination across all shots.

# Requirements
- PlotlyJS.jl must be available for rendering the plots.

# See Also
- [`QuantumOutput`](@ref)

# Example
Plot the measurement outcomes of two qubits entangled in a Bell state:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1000)
julia> probeMeasureOutcome(qo, 4, "Measured outcomes of the 2 measured qubits from a correlated 2-qubit Bell state")
```
"""
function probeMeasureOutcome(quantumOutput::QuantumOutput, stateId::Int, aTitle::String)
    # Ensure stateId is a 1-based index
	stateId = convertToOneBasedNumbering(quantumOutput.settings, stateId)

	# Extract the measurement outcome for each shot
	outcomes = [quantumOutput.output[stateId, k].measured.outcome for k in 1:quantumOutput.numberOfShots]

	# Count the occurrences of each outcome
    measuredLabels = quantumOutput.output[stateId, 1].measured.labels
	orderedOutcomeCounts = [count(outcomes.==(k-1)) for k in 1:length(measuredLabels)] # Assuming outcomes are 0-indexed

	# Normalize the counts to get frequencies
    frequencies = orderedOutcomeCounts / quantumOutput.numberOfShots

	# Create the bar plot
    trace = bar(x=measuredLabels, y=frequencies)
    layout = Layout(
        title_text=aTitle,
        xaxis_title_text="Measured qubit combination",
        yaxis_title_text="Relative Frequency[]"
    )

    # Generate and display the plot
    p = plot([trace], layout)
    display(p)
end

""" 
    probeMeasureProbability(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		aTitle::String
	) -> Nothing

Visualizes the **measured probability distributions** for a specific quantum register in the output of a quantum program.

If the number of shots is 1, this function plots a **bar chart** of the measured probabilities.  
If the number of shots is greater than 1, it plots a **line chart** showing how the probabilities evolve across shots.

# Arguments
- `quantumOutput::QuantumOutput` - The result object returned by `runQuantumProgram`, containing the outcomes of quantum measurements.
- `stateId::Int` - An integer identifying the specific quantum register (subsystem) to probe. It is automatically converted to 1-based indexing if needed.
- `aTitle::String` - The title to display on the histogram plot.

# Returns
- `Nothing` - Makes the plot.

# Behavior
- **Single shot (`numberOfShots == 1`)**:  
  Plots a bar chart of the probabilities associated with each measurement outcome label.
- **Multiple shots (`numberOfShots > 1`)**:  
  Plots a line chart showing how the probability of each possible outcome varies across repeated shots.

# Notes
- `stateId` corresponds to the register index defined during circuit compilation.
- The x-axis will show either the measured label (for bar plots) or shot index (for line plots).
- The y-axis shows the probability (from ``0.0`` to ``1.0``) of each outcome label.

# Requirements
- PlotlyJS.jl must be available for rendering the plots.

# See Also
- [`QuantumOutput`](@ref)

# Example
Plot the measured probabilities for a Bell state:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> measureGate!(qc, [1, 2])
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> probeMeasureProbability(qo, 4, "Measured outcomes of the 2 measured qubits from a correlated 2-qubit Bell state")
```
"""
function probeMeasureProbability(quantumOutput::QuantumOutput, stateId::Int, aTitle::String)         
	# Ensure stateId is a 1-based index
	stateId = convertToOneBasedNumbering(quantumOutput.settings, stateId)

	traces = GenericTrace[]

	if quantumOutput.numberOfShots == 1
		# Extract the measurement outcome for the single shot
		measured = quantumOutput.output[stateId, 1].measured

		# Create the bar plot
		push!(traces, bar(x = measured.labels, y = measured.probability))
		layout = Layout(title_text = aTitle, xaxis_title_text = "Measured qubit combination", yaxis_title_text = "Probability[]")

	else
		# Extract the measurement outcomes for each shot
		measurementsAcrossShots = quantumOutput.output[stateId, :]
		measuredLabels = measurementsAcrossShots[1].measured.labels

		# for each outcome
		for m in 1:length(measuredLabels)

			# Collect the probability of the m-th outcome from each shot
            probabilityOverShots = [shot.measured.probability[m] for shot in measurementsAcrossShots]

			# Create the scatter plot
			push!(traces, scatter(y = probabilityOverShots, mode = "lines", name = measuredLabels[m]))
		end
		layout = Layout(title_text = aTitle, xaxis_title_text = "Shot[#]", yaxis_title_text = "Probability[]")
	end

    # Generate and display the plot
	p = plot(traces, layout)
	display(p)
end

""" 
    probeStateProbability(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		aTitle::String
	) -> Nothing

Visualizes the **state probabilities** of a quantum system (pure or mixed) over one or multiple shots.

This function computes and plots the probability distribution of basis states for a given quantum subsystem specified by `stateId`.
The function handles both `VectorState` (pure states) and `DensityState` (mixed states), adapting the visualization accordingly.

# Arguments
- `quantumOutput::QuantumOutput` - The result of running a quantum circuit, as returned by `runQuantumProgram`.
- `stateId::Int` - The identifier of the quantum register or subsystem to probe. It is automatically converted to 1-based indexing.
- `aTitle::String` - The title to display on the plot(s).

# Returns
- `Nothing` - Makes the plot.

# Behavior
- For **pure states (`VectorState`)**:
  - Computes ``p_m = ‖Mₘ·|ψ⟩‖²`` for each projection operator ``Mₘ``.
- For **mixed states (`DensityState`)**:
  - Computes ``p_m = Tr(Mₘ^†·Mₘ·ρ)`` for each projection operator ``Mₘ``.
- Uses projective measurement operators (Pauli-Z basis) on all qubits to compute the probabilities.

Depending on the number of shots:

- **Single shot (`numberOfShots == 1`)**:  
  A single **bar chart** is displayed showing the probability of each computational basis state.

- **Multiple shots (`numberOfShots > 1`)**:  
  A **line chart** is plotted showing how each state probability evolves over shots.

# Requirements
- This function requires PlotlyJS.jl for rendering.

# See Also
- [`QuantumOutput`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)
- [`calculateEntropyAndPurity`](@ref)

# Example
Plot the state probabilities for the qubits in a Bell state:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> hGate!(qc, 1)
julia> cnotGate!(qc, 1, 2)
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> probeStateProbability(qo, 3, "State probabilities of the 2 qubits from a correlated 2-qubit Bell state")
```
"""
function probeStateProbability(quantumOutput::QuantumOutput, stateId::Int, aTitle::String)
	numberOfQubits = quantumOutput.numberOfQubits
	numberOfShots = quantumOutput.numberOfShots

    # Ensure stateId is a one-based index
    stateId = convertToOneBasedNumbering(quantumOutput.settings, stateId)

    # Generate Kraus operators for PVM measurement
    sigmas = [sigmaZ() for _ in 1:numberOfQubits]
	krausOperators = constructKrausOperatorsForPVMMeasurement(numberOfQubits, Vector(1:numberOfQubits), sigmas)
	labels = Vector([krausOperator.label for krausOperator in krausOperators])

    # Calculate probabilities based on the state type
    probabilities = calculateProbabilities(quantumOutput.output, stateId, krausOperators, numberOfShots)

    # Plot the results
    if numberOfShots == 1
        plotSingleShot(probabilities, labels, aTitle, quantumOutput.output, stateId)
    else
        plotMultiShot(probabilities, labels, aTitle, quantumOutput.output, stateId, numberOfShots)
    end
end

function calculateProbabilities(stateOutput::AbstractMatrix{<:VectorState}, stateId::Int, krausOperators, numberOfShots)
    prob = zeros(length(krausOperators), numberOfShots)

    for (k, shot) in enumerate(stateOutput[stateId,:])
        qVector = shot.q
        for (m, krausOperator) in enumerate(krausOperators)
            Mm = krausOperator.E
            prob[m, k] = sum(abs.(Mm * qVector).^2)
        end
    end
    return prob
end

function calculateProbabilities(stateOutput::AbstractMatrix{<:DensityState}, stateId::Int, krausOperators, numberOfShots)
    probabilities = zeros(length(krausOperators), numberOfShots)

    for (k, shot) in enumerate(stateOutput[stateId,:])
        rhoMatrix = shot.rho
        for (m, krausOperator) in enumerate(krausOperators)
            Mm = krausOperator.E
            probabilities[m, k] = tr((Mm') * Mm * rhoMatrix)
        end
    end
    return probabilities
end

function plotSingleShot(probabilities, labels, aTitle, stateOutput, stateId::Int)
    trace = bar(x=labels, y=probabilities[:, 1])
    layout = Layout(
        title_text=aTitle,
        xaxis_title_text="Qubit combination",
        yaxis_title_text="Probability[]"
    )
    p = plot([trace], layout)
    display(p)
end

function plotMultiShot(probabilities, labels, aTitle, stateOutput, stateId::Int, numberOfShots::Int)
    traces = GenericTrace[]
    for m in 1:size(probabilities, 1)
        push!(traces, scatter(y=probabilities[m, :], mode="lines", name=labels[m]))
    end

    layout = Layout(
        title_text=aTitle,
        xaxis_title_text="Shot[#]",
        yaxis_title_text="Probability[]"
    )
    p = plot(traces, layout)
    display(p)
end

""" 
	probeStateMultiBlochVector(
		quantumOutput::QuantumOutput, 
		stateId::Int, 
		qubits::Qubits, 
		aTitle::String
	) -> Nothing

Visualizes the **Bloch vector representation** of one or more qubits in a quantum state.

This function computes the expectation values ``⟨X⟩``, ``⟨Y⟩``, and ``⟨Z⟩`` for the specified `qubits` from the quantum state identified by `stateId` in `quantumOutput`. The resulting Bloch vectors are plotted on a 3D Bloch sphere.

# Arguments
- `quantumOutput::QuantumOutput` - The result of running a quantum circuit, as returned by `runQuantumProgram`.
- `stateId::Int` - The ID of the subsystem (vector or density state) to be visualized. Internally converted to 1-based indexing.
- `qubits::Qubits` - Indices of the qubits to be plotted, interpreted using the settings for `zeroBasedNumbering` and `indexType`. If empty, all qubits are plotted.
- `aTitle::String` - Title to use for each Bloch sphere visualization.

# Returns
- `Nothing` - Makes the plots.

# Behavior
- Computes the expectation values of Pauli-X, Y, and Z measurements by simulating measurements in each of those bases.
- Displays a translucent Bloch sphere representing the quantum state space for each specified qubit, with:
  - Axes labeled x, y, ``|0⟩``, and ``|1⟩``.
  - The computed Bloch vector shown in green.
  - Coordinate axes shown in red for orientation.

# Requirements
- Requires PlotlyJS.jl for 3D rendering.

# See Also
- [`QuantumOutput`](@ref)
- [`VectorState`](@ref)
- [`DensityState`](@ref)
- [`Qubits`](@ref)

# Example
Plot the Bloch vector for a qubit in a pure superposition state ``|+⟩``:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> hGate!(qc, 1)
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(vector, [1. 0.;])
julia> qo = runQuantumProgram(qp, iqs, 1)
julia> probeStateMultiBlochVector(qo, 2, [1], "Bloch state vector of a qubit in a pure superposition state")
```
"""
function probeStateMultiBlochVector(quantumOutput::QuantumOutput, stateId::Int, qubits::Qubits, aTitle::String)
    numberOfQubits = quantumOutput.numberOfQubits
    settings = quantumOutput.settings

    # --- Data Preparation ---
    stateId = convertToOneBasedNumbering(settings, stateId)
    initialState = quantumOutput.output[stateId, 1] # Assuming single shot analysis
    
    # If no specific qubits are requested, select all of them.
    qubitIndices = isempty(qubits) ? (1:numberOfQubits) : convertToByteIndex(settings, numberOfQubits, qubits)

    # --- Calculation ---
    # This helper function encapsulates the simulation and calculation logic.
    expX, expY, expZ = calculatePauliExpectations(initialState, numberOfQubits)

    # --- Plotting ---
    for qIdx in qubitIndices
        settingsBasedIdx = convertFromOneToSettingsBasedNumbering(settings, qIdx)
        
        println("Qubit $settingsBasedIdx:")
        println("  <X> = $(expX[qIdx])")
        println("  <Y> = $(expY[qIdx])")
        println("  <Z> = $(expZ[qIdx])")

        blochVector = [expX[qIdx], expY[qIdx], expZ[qIdx]]
        plotTitle = "Bloch Vector Qubit $settingsBasedIdx: $aTitle"
        
        # This helper function encapsulates all plotting logic.
        p = createBlochSpherePlot(blochVector, plotTitle)
        display(p)
    end
end

function calculatePauliExpectations(initialState::State, numberOfQubits::Int)
    # Helper to run one type of measurement
    function runMeasurement(pauliOperators)
        qc = createQuantumCircuit(numberOfQubits)
        measureGate!(qc, Vector(1:numberOfQubits), pauliOperators)
        program = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
        # The result state is at index 2 (1=initial, 2=after measurement)
        return runQuantumProgram(program, initialState, 1).output[2, 1].measured
    end

    # Run simulations for all three bases
    measuredX = runMeasurement([sigmaX() for _ in 1:numberOfQubits])
    measuredY = runMeasurement([sigmaY() for _ in 1:numberOfQubits])
    measuredZ = runMeasurement([sigmaZ() for _ in 1:numberOfQubits])

    # Calculate expectation values from probabilities
    expX = zeros(Float64, numberOfQubits)
    expY = zeros(Float64, numberOfQubits)
    expZ = zeros(Float64, numberOfQubits)
    
    numberOfOutcomes = 2^numberOfQubits
    for k in 1:numberOfOutcomes
        # Get the binary representation of the outcome (e.g., [0, 1, 0] for k=3)
        bitValues = reverse(digits(k - 1, base=2, pad=numberOfQubits)')'

        # Eigenvalue is +1 for bit 0, and -1 for bit 1. This is `(-1)^bit_value`.
        eigenvalues = (-1) .^ bitValues

        expX .+= measuredX.probability[k] .* eigenvalues
        expY .+= measuredY.probability[k] .* eigenvalues
        expZ .+= measuredZ.probability[k] .* eigenvalues
    end
    
    return expX, expY, expZ
end

function createBlochSpherePlot(blochVector, title)
    # 1. Create the sphere surface
    n = 40
    u = range(-π, π; length=n)
    v = range(0, π; length=n)
    xSphere = cos.(u) * sin.(v)'
    ySphere = sin.(u) * sin.(v)'
    zSphere = ones(n) * cos.(v)'
    
    sphereSurface = surface(x=xSphere, y=ySphere, z=zSphere,
        showscale=false, opacity=0.3, colorscale="Viridis")

    # 2. Create traces for the axes
    axisLines = [
        # x, y, z axes
        scatter3d(x=[0, 1], y=[0, 0], z=[0, 0], mode="lines", line=attr(width=4, color="grey")),
        scatter3d(x=[0, 0], y=[0, 1], z=[0, 0], mode="lines", line=attr(width=4, color="grey")),
        scatter3d(x=[0, 0], y=[0, 0], z=[0, 1.05], mode="lines", line=attr(width=4, color="grey")),
        scatter3d(x=[0, 0], y=[0, 0], z=[0, -1.05], mode="lines", line=attr(width=4, color="grey")),
    ]

    # 3. Create the Bloch vector trace
    blochVectorTrace = scatter3d(
        x=[0, blochVector[1]], y=[0, blochVector[2]], z=[0, blochVector[3]],
        mode="lines", line=attr(width=8, color="rgb(0, 204, 0)"))

    # 4. Define layout and annotations
    layout = Layout(
        title=title,
        showlegend=false,
        autosize=true,
        scene=attr(
            annotations=[
                attr(showarrow=false, x=1.3, y=0, z=0, text="X", font=attr(size=14)),
                attr(showarrow=false, x=0, y=1.3, z=0, text="Y", font=attr(size=14)),
                attr(showarrow=false, x=0, y=0, z=1.3, text="|0⟩", font=attr(size=16)),
                attr(showarrow=false, x=0, y=0, z=-1.3, text="|1⟩", font=attr(size=16))
            ],
            xaxis=attr(showgrid=true, zeroline=false, showticklabels=false, title=""),
            yaxis=attr(showgrid=true, zeroline=false, showticklabels=false, title=""),
            zaxis=attr(showgrid=true, zeroline=false, showticklabels=false, title="")
        )
    )

    # 5. Combine and return the plot object
    allTraces = [sphereSurface; axisLines; blochVectorTrace]
    return plot(allTraces, layout)
end

end