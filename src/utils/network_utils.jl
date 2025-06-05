module NetworkUtils

include("../constants.jl")
using PyCall
using Graphs
using Statistics
using Random
@pyimport networkx as nx

export load_network, detect_communities, save_network, save_network_analysis_results, save_sensitivity_analysis_results, save_lfr_analysis_results
#-----------------------------------------------------------------------------------------------------------------------------------
function load_network(name::String)
    filename = "$PATH_TO_NETWORKS/$(name).net"
    
    # Verify if the file exists
    if !isfile(filename)
        error("Network does not exist: $filename")
    end
    
    # Load network
    graph = nx.Graph(nx.read_pajek("$PATH_TO_NETWORKS/$name.net"))
    
    # List of communities
    communities = detect_communities(graph)
    
    return graph, communities
end
#-----------------------------------------------------------------------------------------------------------------------------------
function detect_communities(graph)
    communities = nx.algorithms.community.louvain_communities(graph)
    communities = [collect(comm) for comm in communities]
    return communities
end
#-----------------------------------------------------------------------------------------------------------------------------------
function save_network(graph, name::String)
    filename = "$PATH_TO_NETWORKS/$(name).net"
    nx.write_pajek(graph, filename)
end
#-----------------------------------------------------------------------------------------------------------------------------------
function string_keys(dict::Dict{Int, T}) where T
    return Dict(string(k) => v for (k, v) in dict)
end
#-----------------------------------------------------------------------------------------------------------------------------------
end