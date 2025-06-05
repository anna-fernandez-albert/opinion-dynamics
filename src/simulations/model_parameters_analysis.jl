module ModelAnalysis
include("../constants.jl")
include("../utils/network_utils.jl")
include("../model/opinion_dynamics.jl")
include("../visualization/opinion_plots.jl")
using PyCall
using Statistics
using Graphs

using .NetworkUtils
using .OpinionDynamics
using .Visualization

export run_model_sensibility_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
function run_model_sensibility_analysis(graph, communities, name)

    common_plots_path = "$PATH_TO_PLOTS/$(SENSITIVITY_ANALYSIS)"
    size = graph.number_of_nodes()
    NUM_STEPS = 100 * size
    TOLERANCE_STEPS = 5 * size
    
    if isfile("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name).csv")
        println("File already exists, get the previous results")

        results = NetworkUtils.load_network_analysis_results("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name).csv")

    else
        results = []
        open("$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS)_$(name).csv", "a") do file
            write(file, "λ,trust_value,global_period_value,prob_majority_opinion,t_execution,consensus,fraction_positive_final,fraction_negative_final,flip_fraction,constant_positive_fraction,constant_negative_fraction,inestability_dict\n")
        
            for λ_val in λ_VALUES
                for trust_val in LOCAL_TRUST_VALUES
                    for gp_val in GLOBAL_INFLUENCE_PERIODS
                        for pmo_val in PROB_MAJORITY_OPINION_VALUES
                            println("\n\nRunning simulation with λ=$λ_val, trust=$trust_val, gp=$gp_val, pmo=$pmo_val\n")
                            println("The number of steps is: ", NUM_STEPS, " - Number of interactions per node: ", NUM_STEPS / size)

                            initial_opinions = OpinionDynamics.initialize_opinions(communities, OPINION_VALUES, pmo_val)

                            opinion_history, t_execution = OpinionDynamics.run_simulation(graph, initial_opinions, communities, NUM_STEPS, gp_val, λ_val, trust_val, TOLERANCE_STEPS, GLOBAL_INFLUENCE)
                            
                            consensus, fraction_positive_final, fraction_negative_final, flip_fraction, constant_positive_fraction, constant_negative_fraction, inestability_dict = compute_metrics(graph, opinion_history, t_execution, TOLERANCE_STEPS)
                            
                            write(file, "$λ_val,$trust_val,$gp_val,$pmo_val,$t_execution,$consensus,$fraction_positive_final,$fraction_negative_final,$flip_fraction,$constant_positive_fraction,$constant_negative_fraction,$(inestability_dict)\n")

                            push!(results, Dict(
                                "λ" => λ_val,
                                "trust_value" => trust_val,
                                "global_period_value" => gp_val,
                                "prob_majority_opinion" => pmo_val,
                                "t_execution" => t_execution,
                                "consensus" => consensus,
                                "fraction_positive_final" => fraction_positive_final,
                                "fraction_negative_final" => fraction_negative_final,
                                "flip_fraction" => flip_fraction,
                                "constant_positive_fraction" => constant_positive_fraction,
                                "constant_negative_fraction" => constant_negative_fraction,
                                "inestability_dict" => inestability_dict
                            ))
                        end
                    end
                end
            end
        end
    end
    
    # Compute parameter sensibility - correlation between parameters and results - visualization
    correlations = compute_correlations(results)
    println("Correlations: $(correlations)")
    Visualization.plot_correlation_bars(correlations, "$(common_plots_path)/correlation_bars.png")
    Visualization.heatmap_consensus_lambda_tau(results, "$(common_plots_path)/heatmap_consensus_lambda_tau.png")
    Visualization.plot_execution_vs_T(results, "$(common_plots_path)/stability_vs_T.png")
    Visualization.plot_fluctuant_vs_tau(results, "$(common_plots_path)/tau_vs_fluctuant.png")
    Visualization.plot_fluctuant_vs_lambda(results, "$(common_plots_path)/lambda_vs_fluctuant.png")
    Visualization.plot_consensus_vs_tau(results, "$(common_plots_path)/tau_vs_consensus.png")
    Visualization.plot_consensus_vs_lambda(results, "$(common_plots_path)/lambda_vs_consensus.png")

    # Visualize fluctuation analysis for instable nodes
    visualize_fluctuation_analysis(graph, results, common_plots_path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
function compute_metrics(graph, opinion_history, t_execution, tolerance_steps)
    final_state = opinion_history[end, :]

    consensus = abs(mean(final_state))
    fraction_positive_final = mean(final_state .> 0)
    fraction_negative_final = mean(final_state .< 0)

    # For the last tolerance steps, create a dictionary where the key is the node index and the value is a list of the opinion at that step
    flip_analysis = Dict{Int, Vector{Float64}}(node => Float64[] for node in 1:graph.number_of_nodes())
    start_step = max(1, t_execution - tolerance_steps + 1)
    for step in start_step:t_execution
        for node in 1:graph.number_of_nodes()
            push!(flip_analysis[node], opinion_history[step, node])
        end
    end

    # Calculate the fraction of nodes that flipped their opinion in the last tolerance steps
    flip_fraction = 0.0
    constant_positive_fraction = 0.0
    constant_negative_fraction = 0.0
    inestability_dict = Dict{Int, Int}()
    for (node, opinions) in flip_analysis
        all_values_equal = all(opinions .== opinions[1])
        if !all_values_equal
            flip_fraction += 1.0
            inestability_dict[node] = 0
        elseif opinions[1] > 0
            constant_positive_fraction += 1.0
            inestability_dict[node] = 1
        elseif opinions[1] < 0
            constant_negative_fraction += 1.0
            inestability_dict[node] = -1
        end
    end
    println("#Nodes Fliped: $(flip_fraction) - Constant Positive: $(constant_positive_fraction) - Constant Negative: $(constant_negative_fraction) \n")
    flip_fraction /= length(flip_analysis)
    constant_positive_fraction /= length(flip_analysis)
    constant_negative_fraction /= length(flip_analysis)

    # Convert the keys of inestability_dict to strings
    inestability_dict = NetworkUtils.string_keys(inestability_dict)
    return consensus, fraction_positive_final, fraction_negative_final, flip_fraction, constant_positive_fraction, constant_negative_fraction, inestability_dict
end
#-----------------------------------------------------------------------------------------------------------------------------------
function compute_correlations(results)
    λ_vals = [r["λ"] for r in results]
    trust_vals = [r["trust_value"] for r in results]
    global_period_vals = [r["global_period_value"] for r in results]
    pmo_vals = [r["prob_majority_opinion"] for r in results]

    consensus = [r["consensus"] for r in results]

    correlations = Dict(
        "λ" => abs(cor(λ_vals, consensus)),
        "trust_value" => abs(cor(trust_vals, consensus)),
        "global_period_value" => abs(cor(global_period_vals, consensus)),
        "prob_majority_opinion" => abs(cor(pmo_vals, consensus))
    )
    return correlations
end
#-----------------------------------------------------------------------------------------------------------------------------------
function visualize_fluctuation_analysis(graph, results, common_path)
    # Create a directory for the fluctuation analysis results
    mkpath("$common_path/fluctuation_analysis")

    # Filter results for flip_fraction different than 0
    positive_flip_fraction_results = filter(r -> r["flip_fraction"] != 0, results)
    for r in positive_flip_fraction_results
        λ = r["λ"]
        trust_value = r["trust_value"]
        global_period_value = r["global_period_value"]
        prob_majority_opinion = r["prob_majority_opinion"]
        
        analysis = r["inestability_dict"]

        path = "$common_path/fluctuation_analysis/λ_$(λ)_trust_$(trust_value)_gp_$(global_period_value)_pmo_$(prob_majority_opinion).html"

        # Create a graph with each node's color representing the stable opinions (1 for positive, -1 for negative) and the instable opinions (0 for neutral)
        Visualization.plot_network_instable_state(graph, analysis, path)
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
end