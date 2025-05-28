module PolarizationAnalysis

include("../constants.jl")
using PyCall
using Statistics
using LinearAlgebra
using Random
using Distributions
using StatsBase
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.cm as cm
@pyimport seaborn as sns

export calculate_polarization, calculate_opinion_randomnes, 
        analyze_opinion_distribution, track_polarization_over_time
#-----------------------------------------------------------------------------------------------------------------------------------
function calculate_polarization(opinions)
    opinion_values = collect(values(opinions))
    mean_opinion = mean(opinion_values)
    differences = abs.(opinion_values .- mean_opinion)
    polarization = mean(differences)

    max_difference = maximum(differences)
    if max_difference != 0
        polarization /= max_difference
    end
    return polarization
end
#-----------------------------------------------------------------------------------------------------------------------------------
function calculate_opinion_randomnes(opinions)
    opinion_values = collect(values(opinions))

    num_bins = min(10, length(unique(opinion_values)))
    hist = fit(Histogram, opinion_values, nbins=num_bins)
    counts = hist.weights
    
    # Compute proportions and filter out zeros
    proportions = counts ./ sum(counts)
    proportions = filter(p -> p > 0, proportions)
    
    # Shannon entropy
    entropy = -sum(p * log(p) for p in proportions)
    
    max_entropy = log(length(proportions))
    normalized_entropy = entropy / max_entropy
    
    return normalized_entropy
end
#-----------------------------------------------------------------------------------------------------------------------------------
function analyze_opinion_distribution(opinions, name)
    opinion_values = collect(values(opinions))
    
    # Compute polarization metrics
    polarization = calculate_polarization(opinions)
    randomnes = calculate_opinion_randomnes(opinions)
    
    # Visualitzar la distribució d'opinions amb histograma
    plt.figure(figsize=(12, 8))
    
    # Histograma principal
    plt.subplot(2, 2, 1)
    plt.hist(opinion_values, bins=2, color="skyblue", edgecolor="black")
    plt.xlabel("Valors d'Opinió")
    plt.ylabel("Freqüència")
    plt.title("Distribució d'Opinions Finals")
    
    # Afegir mètriques a la visualització
    metrics_text = """
    Polarització: $polarization
    Fragmentació: $randomnes
    """

    plt.subplot(2, 2, 2)
    plt.text(0.1, 0.5, metrics_text, fontsize=12)
    plt.axis("off")
    
    # Gràfic de densitat
    plt.subplot(2, 2, 3)
    sns.kdeplot(opinion_values, fill=true, color="midnightblue")
    plt.xlabel("Valors d'Opinió")
    plt.ylabel("Densitat")
    plt.title("Densitat d'Opinions")
    
    # Boxplot per mostrar la distribució
    plt.subplot(2, 2, 4)
    plt.boxplot(opinion_values, vert=false)
    plt.xlabel("Valors d'Opinió")
    plt.title("Distribució Estadística d'Opinions")
    
    plt.tight_layout()
    plt.savefig("$PATH_TO_PLOTS/opinion_distribution_analysis_$name.png")
    
    # Retornar resultats numèrics
    return Dict(
        "polarization" => polarization,
        "fragmentation" => randomnes
    )
end
#-----------------------------------------------------------------------------------------------------------------------------------
function track_polarization_over_time(opinion_history, name)
    steps = size(opinion_history, 1)
    nodes = size(opinion_history, 2)
    
    polarization_values = zeros(steps)
    randomnes_values = zeros(steps)
    
    # Compute polarization metrics over time
    for t in 1:steps
        temp_opinions = Dict(string(i) => opinion_history[t, i] for i in 1:nodes)
        
        # Compute indexes
        polarization_values[t] = calculate_polarization(temp_opinions)
        randomnes_values[t] = calculate_opinion_randomnes(temp_opinions)
    end
    
    # Visualitzar l'evolució dels índexs
    plt.figure(figsize=(12, 8))
    
    plt.subplot(3, 1, 1)
    plt.plot(1:steps, polarization_values, label="Polarization", color="blue")
    plt.xlabel("Step")
    plt.ylabel("Polarization")
    plt.title("Polarization Index Evolution")
    plt.grid(true)
    
    plt.subplot(3, 1, 2)
    plt.plot(1:steps, randomnes_values, label="Randomness", color="green")
    plt.xlabel("step")
    plt.ylabel("Randomness")
    plt.title("Randomness Index Evolution")
    plt.grid(true)
    
    plt.tight_layout()
    plt.savefig("$PATH_TO_PLOTS/polarization_evolution_$name.png")
    
    return Dict(
        "polarization" => polarization_values,
        "fragmentation" => randomnes_values
    )
end
end