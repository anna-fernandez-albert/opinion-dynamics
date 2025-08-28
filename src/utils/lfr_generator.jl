module LFRGenerator

include("../constants.jl")
using LFRBenchmarkGraphs
using Random
using Graphs
using Statistics

export generate_lfr_networks, load_lfr_network, generate_lfr_network
#-----------------------------------------------------------------------------------------------------------------------------------
function communities_from_cid(cid::Vector{Int})
    # Find unique community IDs, sorted for consistency
    unique_communities = sort(unique(cid))
    
    # Prepare an array of arrays to hold nodes per community
    communities = Vector{Vector{Int}}(undef, length(unique_communities))
    
    for (idx, comm_id) in enumerate(unique_communities)
        # Find all node indices where cid == comm_id
        communities[idx] = findall(x -> x == comm_id, cid)
    end
    
    return communities
end
#-----------------------------------------------------------------------------------------------------------------------------------
function generate_lfr_network(N, μ, k_avg, k_max, name)
    print("Generating LFR network with N=$(N), μ=$(μ), k_avg=$(k_avg), k_max=$(k_max)... ")
    @time lfr, cid = lancichinetti_fortunato_radicchi(N, k_avg, k_max, mixing_parameter=μ)
    lfr_communities = communities_from_cid(cid)

    # Save the network and communities
    save_lfr_network(lfr, lfr_communities, name)

    return lfr, lfr_communities
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
                name = "LFR_N_$(N)_mu_$(μ)_k_avg_$(k_avg)"
                
                # Check if the network already exists
                if isfile("$PATH_TO_NETWORKS/$(name).lgz")
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
function save_lfr_network(lfr, lfr_communities, name, path=PATH_TO_NETWORKS)
    mkpath(path)
    # Save the network in Pajek NET format
    savegraph("$path/$(name).lgz", lfr)

    # Save the communities in a text file
    open("$path/$(name)_communities.txt", "w") do f
        for comm in lfr_communities
            println(f, join(comm, " "))
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function load_lfr_network(name, path)
    # Load network
    lfr = loadgraph("$path/$(name).lgz")

    # Load communities
    lfr_communities = []
    open("$path/$(name)_communities.txt", "r") do f
        for line in eachline(f)
            nodes = parse.(Int, split(chomp(line)))
            push!(lfr_communities, nodes)
        end
    end
    return lfr, lfr_communities
end
end