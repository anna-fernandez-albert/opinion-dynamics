module LFRGenerator

include("../constants.jl")
using LFRBenchmarkGraphs
using PyCall
using Random
using Graphs
using Statistics
@pyimport networkx as nx

export generate_lfr_networks, load_lfr_network
#-----------------------------------------------------------------------------------------------------------------------------------
function generate_lfr_network(N, μ, k_avg, k_max, name)
    name = "lfr_$(name)"
    
    print("Generating LFR network with N=$(N), μ=$(μ), k_avg=$(k_avg), k_max=$(k_max)... ")
    @time lfr, lfr_communities = lancichinetti_fortunato_radicchi(N, k_avg, k_max, mixing_parameter=μ)
    
    graph = nx.Graph()
    
    num_nodes = nv(lfr)
    for i in 1:num_nodes
        graph.add_node(i)
    end
    for edge in edges(lfr)
        graph.add_edge(src(edge), dst(edge))
    end
    
    # Obtain real communities
    communities = []
    for comm_id in unique(lfr_communities)
        comm_nodes = [string(i) for i in 1:num_nodes if lfr_communities[i] == comm_id]
        push!(communities, comm_nodes)
    end

    # Save the network and communities
    save_lfr_network(graph, communities, name)

    return graph, communities
end
#-----------------------------------------------------------------------------------------------------------------------------------
function generate_lfr_networks(N_list, μ_list, k_avg_list)
    networks = Dict()
    
    for N in N_list
        for μ in μ_list
            for k_avg in k_avg_list
                if k_avg < 10 && μ > 0.7
                    println("Skipping generation for N=$(N), μ=$(μ), k_avg=$(k_avg).")
                    continue
                end
                name = "N_$(N)_mu_$(μ)_k_avg_$(k_avg)"
                
                # Check if the network already exists
                if isfile("$PATH_TO_NETWORKS/lfr_$(name).net")
                    println("Network '$name' already exists. Loading it...")
                    graph, communities = load_lfr_network(name, PATH_TO_NETWORKS)
                    networks[(N, μ, k_avg)] = Dict(
                        "network" => graph,
                        "communities" => communities
                    )
                    continue
                end
                
                k_max = min(floor(Int, 3*k_avg), N - 1)
                network, communities = generate_lfr_network(N, μ, k_avg, k_max, name)
                println("Generated network '$name' with N=$(N), μ=$(μ), k_avg=$(k_avg), k_max=$(k_max), and number of communities=$(length(communities))")
                
                networks[(N, μ, k_avg)] = Dict(
                    "network" => network,
                    "communities" => communities
                )
            end
        end
    end
    
    return networks
end
#-----------------------------------------------------------------------------------------------------------------------------------
function save_lfr_network(graph, communities, name, path=PATH_TO_NETWORKS)
    mkpath(path)
    nx.write_pajek(graph, "$path/$(name).net")

    open("$path/$(name)_communities.txt", "w") do file
        for community in communities
            write(file, join(community, " ") * "\n")
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function load_lfr_network(name, path)
    # Load network
    graph = nx.read_pajek("$path/lfr_$(name).net")
    
    if pytypeof(graph) in (nx.MultiGraph, nx.MultiDiGraph)
        println("Converting MultiGraph to simple Graph...")
        graph = nx.Graph(graph)
    end

    # Load communities
    communities = []
    open("$path/lfr_$(name)_communities.txt", "r") do file
        for line in eachline(file)
            push!(communities, split(line))
        end
    end
    
    return graph, communities
end
end