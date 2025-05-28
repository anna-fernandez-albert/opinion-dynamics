module ModelAnalysis
include("../analysis/polarization.jl")
include("../visualization/opinion_plots.jl")
using PyCall
using Statistics
using Graphs

using .PolarizationAnalysis
using .Plots

export run_model_sensitibily_analysis
#-----------------------------------------------------------------------------------------------------------------------------------
const OPINION_VALUE = [-1, 1]
const GLOBAL_INFLUENCE = 1.0
const TRUST_LEVEL = 0.3
const NUM_STEPS = 100000
#-----------------------------------------------------------------------------------------------------------------------------------
function run_model_sensitibily_analysis(graph, communities)
    λ_values = [0.1, 0.3, 0.5, 0.7, 0.9] # Sensitivity to opinion distance influence in the opinion interaction
    trust_values = [0.1, 0.2, 0.3, 0.4, 0.5] # Trust level for local interactions. If the trust level is low, the node will be more influenced by nodes of the same community
    global_period_influence = [500, 1500, 3000, 5000] # Global influence period. The global influence is applied every X steps
    probability_majority_opinion = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]  # Probability of a node adopting the majority opinion of its community
    
    results = Dict()
    for λ_val in λ_values
        for trust_val in trust_values
            for gp_val in global_period_influence
                for pmo_val in probability_majority_opinion
                    PROBABILITY_MAJORITY_OPINION = pmo_val
                    GLOBAL_PERIOD = gp_val
                    TRUST_LEVEL = trust_val
                    λ = λ_val
                    println("Running simulation with λ=$λ_val, trust=$trust_val, gp=$gp_val, pmo=$pmo_val\n\n")

                    opinions = initialize_opinions(communities, PROBABILITY_MAJORITY_OPINION)

                    convergence_step, final_opinions, opinion_history = run_simulation(graph, opinions, communities, NUM_STEPS, GLOBAL_PERIOD, λ, TRUST_LEVEL)
                    
                    # Compute final results
                    final_variance = var(collect(values(final_opinions)))
                    polarization_index = PolarizationAnalysis.calculate_polarization(final_opinions)
                    consensus_level = 1.0 - length(unique(values(final_opinions))) / length(final_opinions)

                    key = (λ_val, trust_val, gp_val, pmo_val)
                    results[key] = Dict(
                        "convergence_step" => convergence_step,
                        "final_variance" => final_variance,
                        "polarization" => polarization_index,
                        "consensus_level" => consensus_level
                    )
                end
            end
        end
    end
    
    Plots.visualize_parameter_sensitivity(results)
    return results
end
#-----------------------------------------------------------------------------------------------------------------------------------
function run_simulation(graph, opinions, communities, max_steps, GLOBAL_PERIOD, λ, TRUST_LEVEL)
    local_opinions = Dict(k => v for (k, v) in opinions)
    
    # Inicialitzar l'historial d'opinions
    size = graph.number_of_nodes()
    opinion_history = zeros(max_steps, size)
    
    convergence_step = max_steps
    
    for step in 1:max_steps
        if different_opinions(local_opinions)
            # Interaccions locals
            node = rand(collect(graph.nodes()))
            local_interaction_update!(graph, local_opinions, node, communities, TRUST_LEVEL, λ)
    
            # Influència global (periòdica)
            if step % GLOBAL_PERIOD == 0
                global_influence_update!(graph, local_opinions, communities, GLOBAL_INFLUENCE, λ, TRUST_LEVEL)
            end
    
            # Registrar opinions
            opinion_values = [local_opinions[string(i)] for i in 1:size]
            opinion_history[step, :] .= opinion_values
        else
            convergence_step = step
            break
        end
    end
    
    return convergence_step, local_opinions, opinion_history
end
#-----------------------------------------------------------------------------------------------------------------------------------
function MVM(communities, opinions, node, λ)
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
function local_interaction_update!(graph, opinions, node, communities, TRUST_LEVEL, λ)
    # Local Interaction based on Multiscale Voter Model
    neighbors = [n for n in graph.neighbors(node)]
    
    neighbor = rand(neighbors)
    Pi = MVM(communities, opinions, node, λ)

    if Pi < TRUST_LEVEL
        opinions[node] = opinions[neighbor]
    end
end
#-----------------------------------------------------------------------------------------------------------------------------------
function global_influence_update!(graph, opinions, communities, GLOBAL_INFLUENCE, λ, TRUST_LEVEL)
    global_opinion = GLOBAL_INFLUENCE

    for node in graph.nodes()
        Pi = MVM(communities, opinions, node, λ)
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
function initialize_opinions(communities, PROBABILITY_MAJORITY_OPINION)
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