module StructureAnalysis

include("../constants.jl")
include("../visualization/opinion_plots.jl")
using PyCall
using Statistics
using LinearAlgebra
using Random
using Graphs
@pyimport networkx as nx

using .Plots

export analyze_network_structure_impact, analyze_bridge_nodes_impact, 
       calculate_structural_metrics, calculate_opinion_structure_correlations,
       identify_bridge_nodes, compare_bridge_regular_dynamics
#-----------------------------------------------------------------------------------------------------------------------------------
function analyze_network_structure_impact(graph, communities, opinion_history, name)
    # Compute structural metrics
    structural_metrics = calculate_structural_metrics(graph, communities)
    
    # Analyze the last step of the opinion dynamics with the structural metrics
    last_step = findlast(row -> length(unique(row)) > 1, eachrow(opinion_history))
    final_opinions = opinion_history[last_step, :]
    
    # Correlation between structural propierties and final opnions
    correlations = calculate_opinion_structure_correlations(structural_metrics, final_opinions)
    
    # Visualize correlations
    Plots.visualize_structure_opinion_correlations(structural_metrics, final_opinions, correlations, communities, name)
    
    return correlations
end
#-----------------------------------------------------------------------------------------------------------------------------------
function calculate_structural_metrics(graph, communities)
    nodes = collect(graph.nodes())
    
    # Inicialitzar diccionari per emmagatzemar les mètriques
    metrics = Dict()
    
    # Centrality metrics
    degree_centrality = nx.degree_centrality(graph)
    betweenness_centrality = nx.betweenness_centrality(graph)
    closeness_centrality = nx.closeness_centrality(graph)

    # Clustering Coefficient
    clustering = nx.clustering(graph)
    
    # Convert to dictionary
    metrics["degree"] = Dict(string(n) => degree_centrality[n] for n in nodes)
    metrics["betweenness"] = Dict(string(n) => betweenness_centrality[n] for n in nodes)
    metrics["closeness"] = Dict(string(n) => closeness_centrality[n] for n in nodes)
    metrics["clustering"] = Dict(string(n) => clustering[n] for n in nodes)
    
    # Compute networks assortativity and modularity
    metrics["assortativity"] = nx.degree_assortativity_coefficient(graph)
    metrics["modularity"] = nx.community.modularity(graph, communities)
    
    return metrics
end
#-----------------------------------------------------------------------------------------------------------------------------------
function calculate_opinion_structure_correlations(structural_metrics, final_opinions)
    n = length(final_opinions)
    
    degree_vec = zeros(n)
    betweenness_vec = zeros(n)
    closeness_vec = zeros(n)
    clustering_vec = zeros(n)
    
    for i in 1:n
        node_str = string(i)
        degree_vec[i] = structural_metrics["degree"][node_str]
        betweenness_vec[i] = structural_metrics["betweenness"][node_str]
        closeness_vec[i] = structural_metrics["closeness"][node_str]
        clustering_vec[i] = structural_metrics["clustering"][node_str]
    end
    
    # Compute correlations
    correlations = Dict()
    correlations["degree"] = cor(degree_vec, final_opinions)
    correlations["betweenness"] = cor(betweenness_vec, final_opinions)
    correlations["closeness"] = cor(closeness_vec, final_opinions)
    correlations["clustering"] = cor(clustering_vec, final_opinions)
    
    return correlations
end
#-----------------------------------------------------------------------------------------------------------------------------------
function analyze_bridge_nodes_impact(graph, communities, opinion_history, name)
    bridge_nodes = identify_bridge_nodes(graph, communities)
    
    compare_bridge_regular_dynamics(bridge_nodes, opinion_history, name)
    
    return bridge_nodes
end
#-----------------------------------------------------------------------------------------------------------------------------------
function identify_bridge_nodes(graph, communities)
    # Dictionary to map nodes to communities
    node_to_community = Dict()
    for (i, comm) in enumerate(communities)
        for node in comm
            node_to_community[node] = i
        end
    end
    
    bridge_nodes = []
    for node in graph.nodes()
        
        neighbor_communities = Set()
        for neighbor in graph.neighbors(node)
            if haskey(node_to_community, neighbor) # Check if neighbor is in the dictionary
                push!(neighbor_communities, node_to_community[neighbor]) # Add the community of the neighbor to the set -> if the community is already in the set, it won't be added again
            end
        end
        
        # If the node is connected to more than one community, it is a bridge node
        if length(neighbor_communities) > 1
            push!(bridge_nodes, node)
        end
    end
    return bridge_nodes
end
#-----------------------------------------------------------------------------------------------------------------------------------
function compare_bridge_regular_dynamics(bridge_nodes, opinion_history, name)
    steps = size(opinion_history, 1)
    nodes = size(opinion_history, 2)
    
    bridge_nodes_int = [parse(Int, node) for node in bridge_nodes]
    bridge_mask = zeros(Bool, nodes)
    bridge_mask[bridge_nodes_int] .= true
    regular_mask = .!bridge_mask
    
    # Average opinions per grup at each time step
    bridge_opinions = [mean(opinion_history[t, bridge_mask]) for t in 1:steps]
    regular_opinions = [mean(opinion_history[t, regular_mask]) for t in 1:steps]
    
    # Variances per grup at each time step
    bridge_variances = [var(opinion_history[t, bridge_mask]) for t in 1:steps]
    regular_variances = [var(opinion_history[t, regular_mask]) for t in 1:steps]
    
    Plots.plot_bridge_regular_dynamics(bridge_opinions, regular_opinions, bridge_variances, regular_variances, name)
end
end