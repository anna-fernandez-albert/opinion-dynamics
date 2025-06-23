module LFRAnalysis

include("../constants.jl")
include("../utils/network_utils.jl")
include("../utils/lfr_generator.jl")
include("../model/opinion_dynamics.jl")
include("../visualization/opinion_plots.jl")

using PyCall
using Statistics
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
    common_plots_path = "$PATH_TO_PLOTS/$(LFR_ANALYSIS)"

    results = []

    # Check if results file exists
    results_file = "$PATH_TO_RESULTS/$(LFR_ANALYSIS).csv"
    existing_combinations = Set{Tuple{Int, Float64, Float64, Float64, Float64}}()

    if isfile(results_file)
        println("File already exists, loading previous results")
        results = NetworkUtils.load_lfr_analysis_results(results_file)
        # Build set of existing parameter combinations
        for r in results
            combo = (r["N"], r["μ"], r["k_avg"], r["λ"], r["trust_level"])
            push!(existing_combinations, combo)
        end
    end

    open(results_file, "a") do file
        # Write header if file is new
        if filesize(results_file) == 0
            write(file, "N,μ,k_avg,modularity,assortativity,clustering_coefficient,std_dev_community_size,av_path_length,λ,trust_level,t_execution,consensus,consensus_per_community,polarization_between_communities,flip_fraction\n")
        end

        for (parameters, network) in networks
            N, μ, k_avg = parameters
            graph, communities = network["network"], network["communities"]

            # Define number of steps and tolerance steps
            NUM_STEPS = 100 * N
            TOLERANCE_STEPS = 5 * N
            GLOBAL_PERIOD = NUM_STEPS / 10

            # Compute Network Properties
            modularity, assortativity, clustering_coefficient, std_dev_community_size, av_path_length = compute_network_properties(graph, communities)

            for λ in λ_VALUES
                for trust_level in LOCAL_TRUST_VALUES
                    combo = (N, μ, k_avg, λ, trust_level)
                    if combo in existing_combinations
                        println("Skipping existing combination: N=$(N), μ=$(μ), k_avg=$(k_avg), λ=$(λ), trust_level=$(trust_level)")
                        continue
                    end

                    println("Running analysis for network with parameters N=$(N), μ=$(μ), k_avg=$(k_avg), λ=$(λ), trust_level=$(trust_level)")

                    # Initialize opinions
                    initial_opinions = OpinionDynamics.initialize_opinions(communities, OPINION_VALUES, PROBABILITY_MAJORITY_OPINION)

                    # Run simulation
                    opinion_history, t_execution = OpinionDynamics.run_simulation(graph, initial_opinions, communities, NUM_STEPS, GLOBAL_PERIOD, λ, trust_level, TOLERANCE_STEPS, GLOBAL_INFLUENCE)

                    consensus, consensus_per_community, polarization_between_communities, flip_fraction = compute_dynamics_results(opinion_history, t_execution, communities, TOLERANCE_STEPS, graph)

                    write(file, "$(N),$(μ),$(k_avg),$(modularity),$(assortativity),$(clustering_coefficient),$(std_dev_community_size),$(av_path_length),$(λ),$(trust_level),$(t_execution),$(consensus),$(consensus_per_community),$(polarization_between_communities),$(flip_fraction)\n")

                    # Store results
                    push!(results, Dict(
                        "N" => N,
                        "μ" => μ,
                        "k_avg" => k_avg,
                        "modularity" => modularity,
                        "assortativity" => assortativity,
                        "clustering_coefficient" => clustering_coefficient,
                        "std_dev_community_size" => std_dev_community_size,
                        "av_path_length" => av_path_length,
                        "λ" => λ,
                        "trust_level" => trust_level,
                        "t_execution" => t_execution,
                        "consensus" => consensus,
                        "consensus_per_community" => consensus_per_community,
                        "polarization_between_communities" => polarization_between_communities,
                        "flip_fraction" => flip_fraction
                    ))
                    # Add to existing combinations to avoid duplicates in the same run
                    push!(existing_combinations, combo)
                end
            end
        end
    end

    correlations = compute_network_correlations(results)
    Visualization.plot_correlation_bars(correlations, "$(common_plots_path)/correlations.png")
    Visualization.plot_consensus_vs_assortativity(results, "$(common_plots_path)/consensus_vs_assortativity.png")
    Visualization.plot_consensus_vs_av_path_length(results, "$(common_plots_path)/consensus_vs_av_path_length.png")
    Visualization.plot_polarization_vs_modularity(results, "$(common_plots_path)/polarization_vs_modularity.png")
    Visualization.plot_consensus_vs_modularity(results, "$(common_plots_path)/consensus_vs_modularity.png")
    Visualization.plot_consensus_vs_std_dev_community_size(results, "$(common_plots_path)/consensus_vs_std_dev_community_size.png")
    Visualization.plot_consensus_vs_clustering_coefficient(results, "$(common_plots_path)/consensus_vs_clustering_coefficient.png")
    Visualization.plot_density_consensus_heatmap(results, "$(common_plots_path)/consensus_heatmap.png")

    # Visualize trust level and λ influence on consensus and stability
    Visualization.plot_consensus_vs_lambda_different_net(results, "$(common_plots_path)/consensus_vs_lambda.png")
    Visualization.plot_consensus_vs_trust_level_different_net(results, "$(common_plots_path)/consensus_vs_trust_level.png")
    Visualization.plot_fluctuant_vs_lambda_different_net(results, "$(common_plots_path)/fluctuant_vs_lambda.png")
    Visualization.plot_fluctuant_vs_trust_level_different_net(results, "$(common_plots_path)/fluctuant_vs_trust_level.png")
    Visualization.heatmap_consensus_lambda_tau(results, "$(common_plots_path)/heatmap_consensus_lambda_tau.png")
    Visualization.plot_consensus_vs_lambda_modularity(results, "$(common_plots_path)/consensus_vs_lambda_modularity.png")
    Visualization.plot_consensus_vs_trust_level_modularity(results, "$(common_plots_path)/consensus_vs_trust_level_modularity.png")
    Visualization.plot_fluctuant_vs_lambda_modularity(results, "$(common_plots_path)/fluctuant_vs_lambda_modularity.png")
    Visualization.plot_fluctuant_vs_trust_level_modularity(results, "$(common_plots_path)/fluctuant_vs_trust_level_modularity.png")
    
end

#-----------------------------------------------------------------------------------------------------------------------------------
function compute_network_properties(graph, communities)
    # Convert Julia communities to Python list of lists with correct node types
    node_type = typeof(first(graph.nodes()))
    if node_type == String
        py_communities = [[String(node) for node in comm] for comm in communities]
    else
        py_communities = [[Int(node) for node in comm] for comm in communities]
    end

    modularity = nx.algorithms.community.quality.modularity(graph, py_communities)
    assortativity = nx.degree_assortativity_coefficient(graph)
    clustering_coefficient = nx.average_clustering(graph)
    std_dev_community_size = std([length(comm) for comm in communities])
    # Calcula la longitud mitjana del camí només si el graf és connex
    if nx.is_connected(graph)
        av_path_length = nx.average_shortest_path_length(graph)
    else
        components = collect(nx.connected_components(graph))
        if isempty(components)
            av_path_length = NaN
        else
            # Troba el conjunt més gran
            largest = components[argmax(map(length, components))]
            subgraph = graph.subgraph(largest)
            av_path_length = nx.average_shortest_path_length(subgraph)
        end
    end

    return modularity, assortativity, clustering_coefficient, std_dev_community_size, av_path_length
end
#-----------------------------------------------------------------------------------------------------------------------------------
function compute_dynamics_results(opinion_history, t_execution, communities, tolerance_steps, graph)
    # Final State of the Simulation
    final_state = opinion_history[end, :]

    # Compute consensus and polarization
    consensus = abs(mean(final_state))
    consensus_per_community = [abs(mean([final_state[parse(Int, node)] for node in comm])) for comm in communities]

    # Polarization between communities
    global_avg = mean(final_state)
    community_means = [mean([final_state[parse(Int, node)] for node in comm]) for comm in communities]
    polarization_between_communities = mean([abs(community_mean - global_avg) for community_mean in community_means])

    # Flip fraction
    # For the last tolerance steps, create a dictionary where the key is the node index and the value is a list of the opinion at that step
    n_nodes = graph.number_of_nodes()
    flip_analysis = Dict{Int, Vector{Float64}}(node => Float64[] for node in 1:n_nodes)
    start_step = max(1, t_execution - tolerance_steps + 1)
    for step in start_step:t_execution
        for node in 1:n_nodes
            push!(flip_analysis[node], opinion_history[step, node])
        end
    end

    # Calculate the fraction of nodes that flipped their opinion in the last tolerance steps
    flip_fraction = 0.0
    for (node, opinions) in flip_analysis
        all_values_equal = all(opinions .== opinions[1])
        if !all_values_equal
            flip_fraction += 1.0
        end
    end
    flip_fraction /= n_nodes
    
    return consensus, consensus_per_community, polarization_between_communities, flip_fraction
end
#-----------------------------------------------------------------------------------------------------------------------------------
function compute_network_correlations(results)
    # Compute correlations between network properties and dynamics results
    correlations = Dict{String, Float64}()

    for key in ["modularity", "assortativity", "clustering_coefficient", "std_dev_community_size", "av_path_length"]
        values = [r[key] for r in results]
        dynamics_values = [r["consensus"] for r in results]
        correlations[key] = cor(values, dynamics_values)
    end

    println("Correlations between network properties and dynamics results:")
    for (key, value) in correlations
        println("  $key: $value")
    end
    return correlations
end
#-----------------------------------------------------------------------------------------------------------------------------------
end