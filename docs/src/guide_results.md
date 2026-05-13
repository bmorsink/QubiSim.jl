# How you process the results

Given a [`QuantumOutput`](@ref), QubiSim offers the following plotting functions to visualize qubit states and measurement outcomes:

1. Histogram of Measured Outcomes: Use [`probeMeasureOutcome`](@ref) to plot the actual measured outcomes for a specific state.
2. Theoretical Measurement Probabilities: Use [`probeMeasureProbability`](@ref) to visualize the theoretical probabilities of measurement outcomes for a specific state.
3. Theoretical Qubit Probabilities: Use [`probeStateProbability`](@ref) to plot the theoretical probabilities of the qubits in a specific state.
4. Bloch Vector Representation: Use [`probeStateMultiBlochVector`](@ref) to represent specific qubits of a state in the Bloch vector format.

Additionally, you can calculate the expectation value of the product of measured qubits — a key metric in correlation experiments — using [`calculateExpectationValueOfProductOfMeasuredQubits`](@ref).
