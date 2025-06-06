module NetworkUtils

include("../constants.jl")
using PyCall
using Graphs
using Statistics
using Random
using DataFrames
using CSV
@pyimport networkx as nx

export load_network, detect_communities, save_network, save_network_analysis_results, save_sensitivity_analysis_results, save_lfr_analysis_results, string_keys, load_network_analysis_results
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
function load_network_analysis_results(filepath)
    # Read CSV with proper type specification
    df = CSV.read(filepath, DataFrame,
                    types=Dict(
                        :λ => Float64,
                        :trust_value => Float64,
                        :global_period_value => Int,
                        :prob_majority_opinion => Float64, 
                        :t_execution => Int,
                        :consensus => Float64,
                        :fraction_positive_final => Float64,
                        :fraction_negative_final => Float64,
                        :flip_fraction => Float64,
                        :constant_positive_fraction => Float64,
                        :constant_negative_fraction => Float64
                    ),
                    silencewarnings=true)
    
    # Convert DataFrame to array of dictionaries
    results = []
    for row in eachrow(df)
        println("Processing row: ", row)
        result = Dict(
            "λ" => row.λ,
            "trust_value" => row.trust_value,
            "global_period_value" => row.global_period_value,
            "prob_majority_opinion" => row.prob_majority_opinion,
            "t_execution" => row.t_execution,
            "consensus" => row.consensus,
            "fraction_positive_final" => row.fraction_positive_final,
            "fraction_negative_final" => row.fraction_negative_final,
            "flip_fraction" => row.flip_fraction,
            "constant_positive_fraction" => row.constant_positive_fraction,
            "constant_negative_fraction" => row.constant_negative_fraction,
            "inestability_dict" => Dict{String, Int}()
        )
        push!(results, result)
    end
    
    println("Successfully loaded $(length(results)) results using CSV.jl")
    return results
end
#-----------------------------------------------------------------------------------------------------------------------------------
function load_lfr_analysis_results(filepath)
    """
    Load LFR analysis results from CSV file using CSV.jl
    """
    # Read CSV with proper type specification
    df = CSV.read(filepath, DataFrame,
                    types=Dict(
                        :N => Int,
                        :μ => Float64,
                        :k_avg => Int,
                        :modularity => Float64,
                        :assortativity => Float64,
                        :clustering_coefficient => Float64,
                        :std_dev_community_size => Float64,
                        :av_path_length => Float64,
                        :λ => Float64,
                        :trust_level => Float64,
                        :t_execution => Int,
                        :consensus => Float64,
                        :consensus_per_community => Float64,
                        :polarization_between_communities => Float64,
                        :flip_fraction => Float64
                    ),
                    silencewarnings=true)
    
    # Convert DataFrame to array of dictionaries
    results = []
    for row in eachrow(df)
        result = Dict(
            "N" => row.N,
            "μ" => row.μ,
            "k_avg" => row.k_avg,
            "modularity" => row.modularity,
            "assortativity" => row.assortativity,
            "clustering_coefficient" => row.clustering_coefficient,
            "std_dev_community_size" => row.std_dev_community_size,
            "av_path_length" => row.av_path_length,
            "λ" => row.λ,
            "trust_level" => row.trust_level,
            "t_execution" => row.t_execution,
            "consensus" => row.consensus,
            "consensus_per_community" => row.consensus_per_community,
            "polarization_between_communities" => row.polarization_between_communities,
            "flip_fraction" => row.flip_fraction
        )
        push!(results, result)
    end
    
    println("Successfully loaded $(length(results)) results using CSV.jl")
    return results

end
#-----------------------------------------------------------------------------------------------------------------------------------
end