# Software architecture and design

## Quantum Circuit

### How QubiSim adds a Unitary gate to a quantum circuit  

- User calls `someUnitaryGate!` with a circuit, qubits and parameters (optional).
- Internally, this:
  - Wraps a factory function, that knows how to create the corresponding unitary operation, into a UnitaryGate.
  - This factory function is typically a closure that captures all the necessary parameters.
  - Computes the byte indices for the target qubits.
  - Calls `addGate!` to insert the gate into the circuit.
    - Decides whether to place the gate in a new `UnitaryStep` (and appends it to the circuit) or merges into an existing `UnitaryStep` (based on bookkeeping).
    - Updates (`numberOfStepsOnQubits`) bookkeeping.
- The factory function later produces the actual unitary operation when needed.

```julia-repl
someUnitaryGate!(quantumCircuit, qubits, p1, ..., pn)
 └── addGate!(quantumCircuit, UnitaryGate(() -> unitaryOperationFactory() /* closure captures parameters p1, ..., pn */, convertToByteIndex(..., [qubits]), "gateName"))
       └── addGate!(…, gate::UnitaryGate)
             ├── based on bookkeeping, determine if new step is needed
             ├── if new step needed → push!(circuit, UnitaryStep([gate]))
             ├── else → push!(existingStep.gates, gate)
             └── update bookkeeping
```

### How QubiSim adds a Measurement gate to a quantum circuit  

- User calls `measureGate!` in one of the following overloads:
  - Just qubits (default Z-basis PVM),
  - Qubits + sigmas (custom measurement basis PVM),
  - Qubits + Kraus operators (arbitrary POVM).
- Internally, this:
  - Creates a `MeasureGate`
  - Calls `addGate!` to insert the gate into the circuit.
    - Places the gate in a new `MeasurementStep` and appends the step to the circuit
    - Updates (`numberOfStepsOnQubits`) bookkeeping.

```julia-repl
measureGate!(quantumCircuit, qubits, sigmas, krausOperators; forgetOutcome=false)
 └── addGate!(quantumCircuit, MeasureGate(convertToByteIndex(..., qubits), sigmas, krausOperators, forgetOutcome))
       └── addGate!(…, gate::MeasureGate)
             ├── push!(circuit, MeasurementStep(gate))
             └── update bookkeeping such that any next new gate insertion invokes the need for a new step
```

### How QubiSim adds a Quantum channel to a quantum circuit  

- User calls `quantumChannelGate!` with qubits + Kraus operators.
- Internally, this:
  - Creates a `QuantumChannelGate`
  - Calls `addGate!` to insert the gate into the circuit.
    - Creates a new `QuantumChannelStep` and appends it to the circuit
    - Updates (`numberOfStepsOnQubits`) bookkeeping.

```julia-repl
quantumChannelGate!(quantumCircuit, qubits, krausOperators)
 ├── barrier!(quantumCircuit)
 ├── addGate!(…, QuantumChannelGate(convertToByteIndex(..., qubits), krausOperators))
 │     └── addGate!(…, gate::QuantumChannelGate)
 │           ├── push!(circuit, QuantumChannelStep(gate))
 │           └── update bookkeeping such that any next new gate insertion invokes the need for a new step
 └── barrier!(quantumCircuit)
```


## Compilation

### How QubiSim compiles a quantum circuit into a quantum program 

- User calls `compileQuantumCircuit` with a `QuantumCircuit`.
- Internally, this:
  - Creates a new empty `QuantumProgram` to hold the compiled representation.
  - Iterates through the circuit, for each step in the circuit, `compileStep!` is invoked.

```julia-repl
compileQuantumCircuit(quantumCircuit; optimizeNumberOfSteps=false)
 ├── QuantumProgram(settings, numberOfQubits, [])
 └── for each step in circuit:
       └── compileStep!(program, step, numberOfQubits, optimize)
```

  - If `compileStep!` is called with a `UnitaryStep`:
    - Starts with an identity `UnitaryOperation`.
    - For each gate:
      - Generates the gate’s unitary matrix using its factory function.
      - Expands it to the full register of qubits.
      - Cascades (multiplies) it into the step’s `UnitaryOperation`.
    - After all gates are processed, either appends or merges this `UnitaryOperation` into the program (depending on optimization).

```julia-repl
compileStep!(…, step::UnitaryStep, …)
 ├── createNQubitOperationId(numberOfQubits)
 ├── for each gate in step.gates:
 │     ├── UGate = gate.unitaryOperationFactory()
 │     ├── moveNQubitOperationToListOfQubits(UGate, gate.qubits, numberOfQubits)
 │     └── cascadeOperations(UStep, UGate)
 ├── if optimize → merge with last UnitaryOperation
 └── else → push!(program.program, UStep)
```

  - If `compileStep!` is called with a `MeasurementStep`:
    - Depending on whether it’s a PVM or POVM measurement:
      - Generates Kraus operators.
    - Expands them to the full register.
    - Depending on forgetOutcome, adds a `MeasureOperation` or `MeasureAndForgetOperation` to the program.

```julia-repl
compileStep!(…, step::MeasurementStep, …)
 ├── if gate.measureGateType == PVM:
 │     ├── sigmaZ() or use gate.sigmas
 │     └── constructKrausOperatorsForPVMMeasurement(numberOfQubits, qubits, sigmas)
 ├── elseif gate.measureGateType == POVM:
 │     └── moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits)
 ├── if gate.forgetOutcome:
 │     └── push!(program.program, MeasureAndForgetOperation(krausOperators))
 └── else:
       └── push!(program.program, MeasureOperation(krausOperators))
```

  - If `compileStep!` is called with a `QuantumChannelStep`:
    - Expands all Kraus operators to the full register.
    - Adds a `QuantumChannelOperation` to the program.

```julia-repl
compileStep!(…, step::QuantumChannelStep, …)
 ├── moveNQubitKrausOperatorToListOfQubits(krausOperator, qubits, numberOfQubits)
 └── push!(program.program, QuantumChannelOperation(krausOperators))
```


## Factory Functions

### Why we use a factory function for unitary gates

- Homogeneous Storage: The `QuantumCircuit` object can hold a simple `Vector{UnitaryGate}`. This is type-stable and efficient. The core logic that iterates through the circuit to simulate it doesn't need a complex if/elseif chain to handle different gates. It just retrieves the unitaryOperationFactory from each `UnitaryGate` and calls it. This is clean and scalable.
- Decoupling: The `QuantumCircuit` and its core processing logic are completely decoupled from the specifics of, e.g., a U3 gate versus a CNOT gate. The circuit only knows about the `UnitaryGate` container, not what's inside it. This makes the system extensible; adding a new gate doesn't require changing the `QuantumCircuit` struct or its simulation loop.
- Lazy Instantiation: The unitary matrix is not created when the gate is added to the circuit. It is created later when the factory function is actually invoked. For large circuits or complex operations, this can be a significant optimization, saving memory and upfront computation.
- Simplicity of the `addGate!` Interface: The `addGate!` function has a single, consistent responsibility. Its job is to place a `UnitaryGate` object into the circuit's data structure. It doesn't need to know how that gate will be realized later.
- Performance: A closure factory is used to capture all necessary gate parameters within the function object itself. This avoids allocating a new data structure for every gate, reducing memory pressure and garbage collection overhead. The Julia compiler specializes and heavily optimizes these closures, resulting in highly performant code at runtime.

```julia-repl
someUnitaryGate!(quantumCircuit, qubits, p1, ..., pn)
 └── addGate!(quantumCircuit, UnitaryGate(() -> unitaryOperationFactory() /* closure captures parameters p1, ..., pn */, convertToByteIndex(..., [qubits]), "gateName"))
     └── (during compilation) gate.unitaryOperationFactory() // here the specialized and optimized version of the closure (for the specific parameters it captured) is called
         └── createNQubitOperationSomeUnitary(p1, ..., pn)
             └── returns unitary matrix for someUnitary operation
```


## Execution

### How QubiSim executes a quantum program on quantum states

- User calls `runQuantumProgram` with:
  - A compiled `QuantumProgram`.
  - An `initialQubitState` (either vector or density representation).
  - A number of shots (repetitions of the program).
- Internally, this:
  - It allocates an output table to hold the quantum state after every step, for every shot.
  - For each step (m, operation) in the program, for each shot (k)
    - Applies the operation to the current quantum state using `applyOperationOnQubitState!`

```julia-repl
runQuantumProgram(quantumProgram, initialQubitState, numberOfShots)
 ├── numberOfStates = length(program) + 1
 ├── output = array of initialQubitState copies
 └── for each (m, operation) in program:
       └── for each shot k:
             └── applyOperationOnQubitState(operation, output[m,k])
                   └── dispatch based on operation + state type
```

  - If `applyOperationOnQubitState!` is called with a `UnitaryOperation` and with a:
    - `VectorState`: multiplies the state vector by the unitary matrix.
    - `DensityState`: performs conjugation 

```julia-repl
applyOperationOnQubitState(unitaryOperation::UnitaryOperation, vectorState::VectorState)
 └── VectorState(unitaryOperation.U * vectorState.q)

applyOperationOnQubitState(unitaryOperation::UnitaryOperation, densityState::DensityState)
 └── DensityState(unitaryOperation.U * densityState.ρ * unitaryOperation.U†)
```

  - If `applyOperationOnQubitState!` is called with a `MeasureOperation` and with a:
    - `VectorState`:
      - Computes probabilities of each measurement outcome from Kraus operators.
      - Randomly selects a possible outcome and collapses state vector to its corresponding pure state.
    - `DensityState`:
      - Computes probabilities of each measurement outcome from Kraus operators.
      - Randomly selects a possible outcome and collapses density matrix to its corresponding pure state.

```julia-repl
applyOperationOnQubitState(measureOperation::MeasureOperation, vectorState::VectorState)
 ├── compute probs from Kraus operators
 ├── randomly select a possible outcome
 └── collapse state vector towards selected outcome and return VectorState with measured info

applyOperationOnQubitState(measureOperation::MeasureOperation, densityState::DensityState)
 ├── compute probs from Kraus operators
 ├── randomly select a possible outcome
 └── collapse density matrix towards selected outcome and return DensityState with measured info
```

  - If `applyOperationOnQubitState!` is called with a `MeasureAndForgetOperation` and with a:
    - `VectorState`: unsupported (error).
    - `DensityState`:
      - Computes probabilities of each measurement outcome from Kraus operators.
      - Creates post-measurement mixed state by summing over all possible outcomes.

```julia-repl
applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, vectorState::VectorState)
 └── throw error (not supported)

applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, densityState::DensityState)
 ├── compute probs from Kraus operators
 └── create post-measurement mixed state by summing over all possible outcomes and return DensityState with nothing measured 
```

  - If `applyOperationOnQubitState!` is called with a `QuantumChannelOperation` and with a:
    - `VectorState`: unsupported (error).
    - `DensityState`: applies Kraus operators and sums over outcomes.

```julia-repl
applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, vectorState::VectorState)
 └── throw error (not supported)

applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, densityState::DensityState)
 ├── apply all Kraus operators and sum over outcomes
 └── return DensityState with nothing measured
```


## Extensibility

### How you can add new functionality to QubiSim

QubiSim uses multiple dispatch in several parts of the codebase, making it straightforward to extend its functionality.

- To add a new unitary gate, we need to define the following functions:
  - `function newUnitaryGate!`
  - `function createNQubitOperationNewUnitaryFactory`
  - `function createNQubitOperationNewUnitary`
- To add a new state, we need to define the following struct and functions:
  - `struct NewState <: State`
  - `function applyOperationOnQubitState(unitaryOperation::UnitaryOperation, newState::NewState)`
  - `function applyOperationOnQubitState(measureOperation::MeasureOperation, newState::NewState)`
  - `function applyOperationOnQubitState(measureAndForgetOperation::MeasureAndForgetOperation, newState::NewState)`
  - `function applyOperationOnQubitState(quantumChannelOperation::QuantumChannelOperation, newState::NewState)`
- To add a new general gate & operation, we need to define the following structs and functions:
  - `struct NewGate <: Gate`
  - `struct NewStep <: Step`
  - `function newGate!`
  - `function addGate!`
  - `function getGate`
  - `function compileStep!`
  - `struct NewOperation <: Operation`
  - `function applyOperationOnQubitState(newOperation::NewOperation, densityState::DensityState)`
  - `function applyOperationOnQubitState(newOperation::NewOperation, vectorState::VectorState)`
