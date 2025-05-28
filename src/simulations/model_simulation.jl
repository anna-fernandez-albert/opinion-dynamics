module OpinionSimulation

include("../constants.jl")
include("../model/opinion_dynamics.jl")
include("../analysis/polarization.jl")
include("../analysis/structure_analysis.jl")
include("../visualization/opinion_plots.jl")

using Random
using Statistics
using .OpinionDynamics
using .PolarizationAnalysis
using .StructureAnalysis
using .Plots

export run_simulation, run_model_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
function run_simulation(graph, communities, opinions, num_steps=NUM_STEPS, global_period=GLOBAL_PERIOD)
    size = graph.number_of_nodes()
    opinion_history = zeros(num_steps, size)
    
    for step in 1:num_steps
        if different_opinions(opinions)
            # Local interactions
            node = rand(collect(graph.nodes()))
            local_interaction_update!(graph, opinions, node, communities)
    
            # Global influence (periodically)
            if step % global_period == 0
                global_influence_update!(graph, opinions, communities)
            end
    
            # Record opinions
            opinion_values = [opinions[string(i)] for i in 1:size]
            opinion_history[step, :] .= opinion_values
        else
            # If all opinions are the same, fill the rest of the history with their value
            for s in step+1:num_steps
                opinion_history[s, :] .= opinion_history[step, :]
            end
            break
        end
    end
    
    return opinion_history
end
#-----------------------------------------------------------------------------------------------------------------------------------
function run_model_analysis(graph, communities, name)
    # Initialize opinions
    opinions = initialize_opinions(communities)
    
    # Run simulation of the opinion dynamics over time
    opinion_history = run_simulation(graph, communities, opinions)
    
    # Analysis of the results
    Plots.plot_opinion_evolution(opinion_history, name)
    Plots.plot_general_variance(opinion_history, name)
    Plots.plot_convergence_time(opinion_history, communities, name)
    Plots.plot_community_trends(opinion_history, communities, name)
    
    # Find the last step where opinions are not all equal
    n_size = graph.number_of_nodes()
    final_step = findlast(row -> length(unique(row)) > 1, eachrow(opinion_history))
    final_opinions = Dict(string(i) => opinion_history[final_step, i] for i in 1:n_size)
    
    # Compare with inicial step and intermediate step
    # Perform analysis
    polarization_index = PolarizationAnalysis.analyze_opinion_distribution(final_opinions, name)
    polarization_over_time = PolarizationAnalysis.track_polarization_over_time(opinion_history, name)
    structure_correlations = StructureAnalysis.analyze_network_structure_impact(graph, communities, opinion_history, name)
    bridges = StructureAnalysis.analyze_bridge_nodes_impact(graph, communities, opinion_history, name)
    
    # Retornar resultats en un diccionari
    return Dict(
        "opinion_history" => opinion_history,
        "final_opinions" => final_opinions,
        "polarization_index" => polarization_index,
        "polarization_over_time" => polarization_over_time,
        "structure_correlations" => structure_correlations,
        "bridge_nodes" => bridges,
        "convergence_step" => final_step
    )
end

end