module OpinionDynamics

using Statistics
using Distributions
using Random
using Graphs

export run_simulation, MVM, local_interaction_update!, global_influence_update!, different_opinions, initialize_opinions, initialize_trust_levels
#-----------------------------------------------------------------------------------------------------------------------------------
function MVM(communities, opinions, node, λ)
    Pi = 0
    community_opinion = 0.0

    # Compute community opinion (s_c)
    for comm in communities
        if node in comm
            if length(comm) > 1
                community_opinion = sum([opinions[l] for l in comm if l != node]) / (length(comm) - 1)
            else
                community_opinion = opinions[node]
            end
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
function local_interaction_update!(graph, opinions, node, communities, λ, trust_level)
    # Local Interaction based on Multiscale Voter Model
    node_neighbors = neighbors(graph, node)
    if isempty(node_neighbors)
        return  # No neighbors to interact with
    end

    neighbor = rand(node_neighbors)
    Pi = MVM(communities, opinions, node, λ)

    if rand() < Pi * trust_level
        opinions[node] = opinions[neighbor]
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function global_influence_update!(graph, opinions, communities, λ, trust_level, global_opinion)
    # Global Influence based on Multiscale Voter Model
    for node in graph.nodes()
        Pi = MVM(communities, opinions, node, λ)
        if rand() < Pi * trust_level
            opinions[node] = global_opinion
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function different_opinions(opinion_history, step, tolerance_steps; epsilon=1e-2)
    if step <= tolerance_steps
        return true
    end
    recent = opinion_history[(step - tolerance_steps + 1):step, :]
    base = recent[1, :]
    for i in range(2, size(recent, 1))
        if any(abs.(recent[i, :] .- base) .> epsilon)
            return true # There are different opinions in the recent steps
        end
    end
    return false # All recent steps have the same opinion
end
#-----------------------------------------------------------------------------------------------------------------------------------
function initialize_opinions(size, communities, opinion_values, prob_majority_opinion)
    # Array of opinions (1-based, integer indices): opinions[i] gives the opinion of node i
    opinions = zeros(Int, size)
    for community in communities
        majority_opinion = rand(opinion_values)
        for node in community
            if rand() < prob_majority_opinion
                # Assign the majority opinion to a certain percentage of nodes
                opinions[node] = majority_opinion
            else
                opinions[node] = -majority_opinion  # Assign the opposite opinion to the rest
            end
        end
    end
    return opinions
end
#-----------------------------------------------------------------------------------------------------------------------------------
function run_simulation(graph, opinions, communities, num_steps, λ, trust_level, tolerance_steps)    
    """
    Run the simulation of opinion dynamics on the given graph.
    Returns the opinion history and the number of steps taken.

    Parameters:
        - graph: NetworkX graph object representing the social network.
        - opinions: Array of initial opinions for each node. 1D array of integers.
        - communities: List of communities, where each community is a list of node indices.
        - num_steps: Number of steps to run the simulation.
        - global_period: Period for global influence updates.
        - λ: Parameter for the Multiscale Voter Model.
        - trust_level: Trust level for interactions.
        - tolerance_steps: Number of steps to consider for convergence.
        - global_influence_opinion: Opinion to be adopted during global influence updates.
    """
    size = nv(graph)
    #opinions = deepcopy(initial_opinions) # Unique vector of opinions to update
    opinion_history = zeros(num_steps+1, length(communities)) # num_steps x num communities matrix to store the history of opinions
    opinion_history[1, :] .= [(sum(opinions[i] == 1 for i in com) / length(com)) for com in communities]
    println("Opinion history initialized, initial opinions: ", opinion_history[1, :])
    
    for step in 2:num_steps+1 #TODO: epochs = 1000 -> nodes permutation per epoch. -> at each epoch, random order of nodes (permutation) and each node updates its opinion based on the neighbors.
        println("Step: $step")

        # Randomly permute the order of nodes for interaction
        nodes_order = randperm(size)
        for node in nodes_order
            # Local interactions
            OpinionDynamics.local_interaction_update!(graph, opinions, node, communities, λ, trust_level)
        end
        # Save variables to analyze: for each community, save the number of nodes that have opinion 1.
        opinion_history[step, :] .= [(sum(opinions[i] == 1 for i in com) / length(com)) for com in communities]

        # Check if the opinions have converged
        if !OpinionDynamics.different_opinions(opinion_history, step, tolerance_steps)
            return opinion_history[1:step,:], step
        end
    end
    return opinion_history, num_steps
end
#-----------------------------------------------------------------------------------------------------------------------------------
end