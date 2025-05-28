module OpinionDynamics

include("../constants.jl")
using Statistics
using Random

export MVM, local_interaction_update!, global_influence_update!, different_opinions, initialize_opinions, initialize_trust_levels
#-----------------------------------------------------------------------------------------------------------------------------------
function MVM(communities, opinions, node)
    Pi = 0
    community_opinion = 0.0

    # Compute community opinion
    for comm in communities
        if node in comm
            community_opinion = sum([opinions[l] for l in comm if l != node]) / (length(comm) - 1)
            break
        end
    end

    # Compute the opinion distance between the node and the community
    opinion_distance = abs(opinions[node] - community_opinion)

    # Probability of node i adopting the state of a neighbor
    Pi = 1 / (1 + exp((1 - opinion_distance) / λ))
    return Pi
end
#-----------------------------------------------------------------------------------------------------------------------------------
function local_interaction_update!(graph, opinions, node, communities)
    # Local Interaction based on Multiscale Voter Model
    neighbors = [n for n in graph.neighbors(node)]
    
    neighbor = rand(neighbors)
    Pi = MVM(communities, opinions, node)

    if Pi < TRUST_LEVEL
        opinions[node] = opinions[neighbor]
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function global_influence_update!(graph, opinions, communities)
    global_opinion = GLOBAL_INFLUENCE

    for node in graph.nodes()
        Pi = MVM(communities, opinions, node)
        if Pi < TRUST_LEVEL
            opinions[node] = global_opinion
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function different_opinions(opinions)
    return length(unique(values(opinions))) > 1
end
#-----------------------------------------------------------------------------------------------------------------------------------
function initialize_opinions(communities)
    opinions = Dict{String, Float64}()
    for community in communities
        majority_opinion = rand(OPINION_VALUE)
        for node in community
            if rand() < PROBABILITY_MAJORITY_OPINION
                # Assign the majority opinion to a certain percentage of nodes
                opinions[node] = majority_opinion
            else
                opinions[node] = rand(OPINION_VALUE)
            end
        end
    end
    return opinions
end
#-----------------------------------------------------------------------------------------------------------------------------------
function initialize_trust_levels(graph)
    positions = nx.spring_layout(graph)
    maximum_distance = maximum([norm(positions[i] - positions[j]) for i in graph.nodes(), j in graph.nodes()]) # or radious

    # Pre-allocate trust levels matrix
    trust_levels = Dict(node => Dict() for node in graph.nodes())

    for node in graph.nodes()
        neighbors = [n for n in graph.neighbors(node)]
        for neighbor in neighbors
            normalized_distance = norm(positions[node] - positions[neighbor]) / maximum_distance
            trust_levels[node][neighbor] = normalized_distance
        end
    end   
    return trust_levels
end
#-----------------------------------------------------------------------------------------------------------------------------------
end