module LFRAnalysis

include("../constants.jl")
include("../utils/network_utils.jl")
include("../utils/lfr_generator.jl")
include("../model/opinion_dynamics.jl")
include("../visualization/opinion_plots.jl")

using PyCall
using Statistics
using Graphs
using Graphs: modularity, assortativity, global_clustering_coefficient
using CSV, DataFrames
@pyimport networkx as nx

using .Visualization
using .NetworkUtils
using .LFRGenerator
using .OpinionDynamics

export run_lfr_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
function run_lfr_analysis()
    # Create networks
    networks = LFRGenerator.generate_lfr_networks(N_VALUES, μ_VALUES, k_avg_VALUES)

    for (parameters, network) in networks
        N, μ, k_avg = parameters
        name = "LFR_N_$(N)_mu_$(μ)_k_avg_$(k_avg)"
        graph, communities = network["network"], network["communities"]
        println("\n\nAnalyzing network with parameters N=$(N), μ=$(μ), k_avg=$(k_avg) - Name: $name\n")

        # Define number of steps and tolerance steps
        size = nv(graph)
        repetitions = size < 5000 ? 100 : 20

        if isfile("$PATH_TO_NETWORKS/$(name)_properties.txt")
            println("File already exists - we do not compute network properties again")
        else
            # Compute Network Properties
            modularity, assortativity, clustering_coefficient, std_dev_community_size = compute_network_properties(graph, communities)
            
            # Save to TXT file
            open("$PATH_TO_NETWORKS/$(name)_properties.txt", "w") do file
                write(file, "modularity,assortativity,clustering_coefficient,std_dev_community_size\n")
                write(file, "$(modularity),$(assortativity),$(clustering_coefficient),$(std_dev_community_size)\n")
            end
        end

        for λ in λ_VALUES
            for trust_level in LOCAL_TRUST_VALUES
                println("\n\nRunning simulation with trust=$trust_level, λ=$λ\n")
                for i in 1:repetitions
                    println("Repetition $i of $repetitions")
                    print("Parameters: N=$(N), μ=$(μ), k_avg=$(k_avg) -- trust=$(trust_level), λ=$(λ)\n")
                    if isfile("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name)_parameters_$(PROBABILITY_MAJORITY_OPINION)_$(trust_level)_$(λ)_repetition_$(i).csv")
                        println("File already exists")
                    else
                        current_time = time()
                        initial_opinions = OpinionDynamics.initialize_opinions(size, communities, OPINION_VALUES, PROBABILITY_MAJORITY_OPINION)

                        opinion_history, t_execution = OpinionDynamics.run_simulation(graph, initial_opinions, communities, NUM_STEPS, λ, trust_level, TOLERANCE_STEPS)
                        print("Execution time: $t_execution steps")
                        df = DataFrame(opinion_history, :auto)
                        CSV.write("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name)_parameters_$(PROBABILITY_MAJORITY_OPINION)_$(trust_level)_$(λ)_repetition_$(i).csv", df)
                        print("Time taken: $(time() - current_time) seconds")
                    end
                end
            end
        end
    end
end

#-----------------------------------------------------------------------------------------------------------------------------------
function compute_network_properties(graph, comms)
    # Modularity
    # Modularity (Graphs.jl expects a partition vector)
    partition = zeros(Int, nv(graph))
    for (i, comm) in enumerate(comms)
        for node in comm
            partition[node] = i
        end
    end
    mod = modularity(graph, partition)

    # Assortativity
    assort = assortativity(graph)

    # Clustering coefficient (average over all nodes)
    clust_coeff = mean(global_clustering_coefficient(graph))

    # Standard deviation of community sizes
    std_dev_comm_size = std([length(comm) for comm in comms])

    return mod, assort, clust_coeff, std_dev_comm_size
end
#-----------------------------------------------------------------------------------------------------------------------------------
end