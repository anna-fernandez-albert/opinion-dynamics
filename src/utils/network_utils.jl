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
function save_network_analysis_results(results, name::String)
    filename = "$PATH_TO_RESULTS/$(NETWORK_ANALYSIS_RESULTS_FILE_NAME).csv"
    
    # if exists do not write header
    if isfile(filename)
        open(filename, "a") do file
            write(file, "$name,$(results["polarization_index"]),$(results["structure_correlations"]["degree"]),$(results["structure_correlations"]["betweenness"]),$(results["structure_correlations"]["closeness"]),$(results["structure_correlations"]["clustering"]),$(results["bridge_nodes"]),$(results["convergence_step"])")
        end
    else
        open(filename, "w") do file
            write(file, "network_name,polarization_index,degree_correlation,betweenness_correlation,closeness_correlation,clustering_correlation,bridge_nodes,convergence_step\n")
            write(file, "$name,$(results["polarization_index"]),$(results["structure_correlations"]["degree"]),$(results["structure_correlations"]["betweenness"]),$(results["structure_correlations"]["closeness"]),$(results["structure_correlations"]["clustering"]),$(results["bridge_nodes"]),$(results["convergence_step"])")
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function save_sensitivity_analysis_results(results, name::String)
    filename = "$PATH_TO_RESULTS/$(SENSITIVITY_ANALYSIS_RESULTS_FILE_NAME).csv"

    open(filename, "a") do file
        # Write header
        write(file, "λ,trust_value,global_period_value,prob_majority_opinon,convergence_step,final_variance,polarization,consensus_level\n")
        
        for key in keys(results)
            λ_val, trust_val, gp_val, gm_val = key
            write(file, "$λ_val,$trust_val,$gp_val,$gm_val,$(results[key]["convergence_step"]),$(results[key]["final_variance"]),$(results[key]["polarization"]),$(results[key]["consensus_level"])\n")    
        end
    end
    
end
#-----------------------------------------------------------------------------------------------------------------------------------
function save_lfr_analysis_results(results, name::String)
    filename = "$PATH_TO_RESULTS/$(name)_lfr_analysis.csv"

    open(filename, "a") do file
        # Write header
        write(file, "N,μ,k_avg,polarization_index,degree_correlation,betweenness_correlation,closeness_correlation,clustering_correlation,bridge_nodes,convergence_step\n")
        
        for key in keys(results)
            N_val, μ_val, k_avg_val = key
            write(file, "$N_val,$μ_val,$k_avg_val,$(results[key]["polarization_index"]),$(results[key]["structure_correlations"]["degree"]),$(results[key]["structure_correlations"]["betweenness"]),$(results[key]["structure_correlations"]["closeness"]),$(results[key]["structure_correlations"]["clustering"]),$(results[key]["bridge_nodes"]),$(results[key]["convergence_step"])")
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
end