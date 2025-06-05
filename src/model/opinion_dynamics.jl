module OpinionDynamics

using Statistics
using Distributions
using Random

export run_simulation, MVM, local_interaction_update!, global_influence_update!, different_opinions, initialize_opinions, initialize_trust_levels
#-----------------------------------------------------------------------------------------------------------------------------------
function MVM(communities, opinions, node, λ)
    Pi = 0
    community_opinion = 0.0

    # Compute community opinion (s_c)
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
function local_interaction_update!(graph, opinions, node, communities, λ, trust_level)
    # Local Interaction based on Multiscale Voter Model
    neighbors = collect(graph.neighbors(node))
    
    neighbor = rand(neighbors)
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
function different_opinions(opinion_history, step, tolerance_steps)
    if step <= tolerance_steps
        return true
    end
    recent = opinion_history[(step - tolerance_steps + 1):step, :]
    # Comprova si totes les files són iguals entre si
    for i in 2:size(recent, 1)
        if any(recent[i, :] .!= recent[1, :])
            return true  # Hi ha hagut canvi d'opinió
        end
    end
    return false  # No hi ha hagut canvi d'opinió
end
#-----------------------------------------------------------------------------------------------------------------------------------
function initialize_opinions(communities, opinion_values, prob_majority_opinion)
    opinions = Dict{String, Float64}()
    for community in communities
        majority_opinion = rand(opinion_values)
        for node in community
            if rand() < prob_majority_opinion
                # Assign the majority opinion to a certain percentage of nodes
                opinions[node] = majority_opinion
            else
                opinions[node] = rand(opinion_values)
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
function run_simulation(graph, initial_opinions, communities, num_steps, global_period, λ, trust_level, tolerance_steps, global_influence_opinion)    
    size = graph.number_of_nodes()
    opinions = deepcopy(initial_opinions)
    opinion_history = zeros(num_steps, size)
    
    for step in 1:num_steps
        # Record opinions
        opinion_history[step, :] .= [opinions[string(i)] for i in 1:size]

        if OpinionDynamics.different_opinions(opinion_history, step, tolerance_steps)
            # Local interactions
            node = rand(collect(graph.nodes()))
            OpinionDynamics.local_interaction_update!(graph, opinions, node, communities, λ, trust_level)
    
            # Global influence (periodically)
            if step % global_period == 0
                OpinionDynamics.global_influence_update!(graph, opinions, communities, λ, trust_level, global_influence_opinion)
            end
        else
            # If opinions are converged, fill the rest of the history with the last opinion
            for s in step+1:num_steps
                opinion_history[s, :] .= opinion_history[step, :]
            end

            return opinion_history, step
        end
    end
    
    return opinion_history, num_steps
end
#-----------------------------------------------------------------------------------------------------------------------------------
end