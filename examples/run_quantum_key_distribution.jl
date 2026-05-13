using QubiSim
using LinearAlgebra, PlotlyJS
using Random

function encode_message(bits, bases)
    message = []
    for i in 1:n
        qc=createQuantumCircuit(1)
        if bases[i]==0 # Prepare qubit in Z-basis
            if bits[i]==1
                xGate!(qc, 1) 
            else
                idGate!(qc, 1)
            end 
        else # Prepare qubit in X-basis
            if bits[i]==1
                xGate!(qc, 1) 
                hGate!(qc, 1) 
            else
                hGate!(qc, 1)
            end
        end
        barrier!(qc)
        push!(message,qc)
    end
    return message
end

function measure_message!(message, bases)
    measurements = []
    for q in 1:n
        if bases[q]==0 # measuring in Z-basis
            measureGate!(message[q], [1])
        end
        if bases[q]==1 # measuring in X-basis
            hGate!(message[q], 1)
            measureGate!(message[q], [1])
        end

        qp=compileQuantumCircuit(message[q]; optimizeNumberOfSteps=true)
        iqs=createInitialQubitState(vector, [1. 0.])
        qo=runQuantumProgram(qp, iqs, 1)

        measured_bit=qo.output[end].measured.outcome[1] # note [end] is important otherwise we will sample eve's interception instead..
        push!(measurements, measured_bit)
    end
    return measurements
end

function remove_garbage(a_bases, b_bases, bits)
    good_bits = []
    for q in 1:n
        if a_bases[q] == b_bases[q]
            # If both used the same basis, add this to the list of 'good' bits
            push!(good_bits, bits[q])
        end
    end
    return good_bits
end

function sample_bits!(bits, selection)
    sample = []
    for i in selection
        # use np.mod to make sure the bit we sample is always in the list range
        i = mod(i, length(bits))+1
        bit=bits[i]
        deleteat!(bits,i)
        push!(sample, bit)
    end
    return sample
end

print("\n")
print("1. Quantum key distribution: without interception\n")

n = 100
## Step 1
# Alice generates bits
alice_bits = rand(0:1, n)

## Step 2
# Create an array to tell us which qubits are encoded in which bases
alice_bases = rand(0:1, n)
message = encode_message(alice_bits, alice_bases)

## Step 3
# Decide which basis to measure in:
bob_bases = rand(0:1, n)
bob_results = measure_message!(message, bob_bases)

## Step 4
alice_key = remove_garbage(alice_bases, bob_bases, alice_bits)
bob_key = remove_garbage(alice_bases, bob_bases, bob_results)

## Step 5
sample_size = 15
bit_selection = rand(1:n, sample_size)
bob_sample = sample_bits!(bob_key, bit_selection)
alice_sample = sample_bits!(alice_key, bit_selection)
print("  bob_sample = ", bob_sample',"\n")
print("alice_sample = ", alice_sample',"\n")
print("bob_sample==alice_sample : ",bob_sample==alice_sample,"\n")
print("  bob_key = ", bob_key',"\n")
print("alice_key = ", alice_key',"\n")
print("key length = ", length(alice_key),"\n")


print("\n")
print("2. Quantum key distribution: with interception\n")

n = 200
## Step 1
# Alice generates bits
alice_bits = rand(0:1, n)

## Step 2
# Create an array to tell us which qubits are encoded in which bases
alice_bases = rand(0:1, n)
message = encode_message(alice_bits, alice_bases)

## Interception!!
eve_bases = rand(0:1, n)
intercepted_message = measure_message!(message, eve_bases)
# be careful, in this way eve will get the intercepted message but later when bob measures, the circuit is extended with
# bob's measurement apparatus and the whole circuit is rerun, so another intercepted message will be retrieved.

## Step 3
# Decide which basis to measure in:
bob_bases = rand(0:1, n)
bob_results = measure_message!(message, bob_bases)

## Step 4
alice_key = remove_garbage(alice_bases, bob_bases, alice_bits)
bob_key = remove_garbage(alice_bases, bob_bases, bob_results)

print("probability bits equal in both keys=",sum(alice_key.==bob_key)/length(alice_key),"\n")

## Step 5
sample_size = 15
bit_selection = rand(1:n, sample_size)
bob_sample = sample_bits!(bob_key, bit_selection)
alice_sample = sample_bits!(alice_key, bit_selection)
print("  bob_sample = ", bob_sample',"\n")
print("alice_sample = ", alice_sample',"\n")
print("bob_sample==alice_sample : ",bob_sample==alice_sample,"\n")
print("  bob_key = ", bob_key',"\n")
print("alice_key = ", alice_key',"\n")
print("key length = ", length(alice_key),"\n")
