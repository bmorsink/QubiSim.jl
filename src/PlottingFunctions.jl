############################################
# PlottingFunctions

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
