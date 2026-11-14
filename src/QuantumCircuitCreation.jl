############################################
# QuantumCircuitCreation

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
function addGate!(quantumCircuit::QuantumCircuit, gate::BernoulliGate)
	push!(quantumCircuit.circuit, BernoulliStep(gate))
	quantumCircuit.numberOfStepsOnQubits = (length(quantumCircuit.circuit)+1)*ones(Int, quantumCircuit.numberOfQubits)
end
function addGate!(quantumCircuit::QuantumCircuit, gate::StochasticGate)
	push!(quantumCircuit.circuit, StochasticStep(gate))
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