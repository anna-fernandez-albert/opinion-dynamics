module Plots

include("../constants.jl")
using PyCall
using Statistics
using LinearAlgebra
@pyimport networkx as nx
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.cm as cm
@pyimport pyvis.network as pyvis

export plot_network, plot_opinion_evolution, plot_community_trends, plot_general_variance, 
       plot_convergence_time, visualize_structure_opinion_correlations, plot_bridge_regular_dynamics,
       visualize_parameter_sensitivity
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_network(graph, communities, opinions, name)
    G = pyvis.Network(height="1000px", width="100%")
    pos = nx.spring_layout(graph)
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray", "cyan", "magenta", "yellow", "black", "white"]
    j = 1
    for comm in communities
        for i in comm
            G.add_node(i, x=pos[i][1]*10000, y=pos[i][2]*10000, size=150, physics=false, color=colors[j], label=opinions[i])
        end
        j += 1
    end
    for i in graph.edges
        G.add_edge(i[1], i[2], physics=false, width=0.5, color="gray")     
    end
    G.save_graph("$PATH_TO_PLOTS/$name.html")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_opinion_evolution(opinion_history, name)
    steps = size(opinion_history, 1)
    avg_opinions = [mean(opinion_history[t, :]) for t in 1:steps]

    plt.figure(figsize=(10, 6))
    plt.plot(1:steps, avg_opinions, label="Network Average Opinion", color="midnightblue")
    plt.xlabel("Time Step")
    plt.ylabel("Average Opinion")
    plt.legend()
    plt.savefig("$PATH_TO_PLOTS/opinion_evolution_$name.png")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_community_trends(opinion_history, communities, name)
    steps = size(opinion_history, 1)
    community_avg = [mean([opinion_history[t, parse(Int, node)] for node in comm]) for comm in communities, t in 1:steps]

    plt.figure(figsize=(10, 6))
    for i in 1:length(communities)
        plt.plot(1:steps, community_avg[i, :], label="Community $i", linestyle="--", color=cm.Blues(i/ length(communities)))
    end
    plt.xlabel("Time Step")
    plt.ylabel("Community Average Opinion")
    plt.legend()
    plt.savefig("$PATH_TO_PLOTS/community_trends_$name.png")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_general_variance(opinion_history, name)
    steps = size(opinion_history, 1)
    general_variances = [var(opinion_history[t, :]) for t in 1:steps]

    plt.figure(figsize=(10, 6))
    plt.plot(1:steps, general_variances, label="Network Variance", color="midnightblue")
    plt.xlabel("Time Step")
    plt.ylabel("Network Variance")
    plt.legend()
    plt.savefig("$PATH_TO_PLOTS/general_variance_$name.png")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_convergence_time(opinion_history, communities, name)
    steps = size(opinion_history, 1)
    community_variances = []

    for comm in communities
        variances = []
        for t in 1:steps
            opinions_comm = [opinion_history[t, parse(Int, node)] for node in comm]
            push!(variances, var(opinions_comm))
        end
        push!(community_variances, variances)
    end

    # Plot community variances over time
    plt.figure(figsize=(10, 6))
    for i in 1:length(communities)
        plt.plot(1:steps, community_variances[i], label="Community $i", linestyle="--", color=cm.Blues(i/ length(communities)))
    end
    plt.xlabel("Time Step")
    plt.ylabel("Community Variance")
    plt.legend()
    plt.savefig("$PATH_TO_PLOTS/community_variances_$name.png")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function visualize_structure_opinion_correlations(structural_metrics, final_opinions, correlations, communities, name)
    """
    Visualitza les correlacions entre propietats estructurals i opinions finals.
    """
    n = length(final_opinions)
    
    # Crear vectors per a les propietats estructurals
    degree_vec = zeros(n)
    betweenness_vec = zeros(n)
    closeness_vec = zeros(n)
    clustering_vec = zeros(n)
    
    for i in 1:n
        node_str = string(i)
        degree_vec[i] = structural_metrics["degree"][node_str]
        betweenness_vec[i] = structural_metrics["betweenness"][node_str]
        closeness_vec[i] = structural_metrics["closeness"][node_str]
        clustering_vec[i] = structural_metrics["clustering"][node_str]
    end
    
    # Crear subplots per a cada propietat estructural
    plt.figure(figsize=(15, 12))
    
    # Degree centrality vs opinion
    plt.subplot(2, 3, 1)
    plt.scatter(degree_vec, final_opinions, alpha=0.5)
    plt.xlabel("Degree Centrality")
    plt.ylabel("Final Opinion")
    plt.title("Degree vs Opinion (r = $(correlations["degree"]))")
    plt.grid(true)
    
    # Betweenness centrality vs opinion
    plt.subplot(2, 3, 2)
    plt.scatter(betweenness_vec, final_opinions, alpha=0.5, color="orange")
    plt.xlabel("Betweenness Centrality")
    plt.ylabel("Final Opinion")
    plt.title("Betweenness vs Opinion (r = $(correlations["betweenness"]))")
    plt.grid(true)
    
    # Closeness centrality vs opinion
    plt.subplot(2, 3, 3)
    plt.scatter(closeness_vec, final_opinions, alpha=0.5, color="green")
    plt.xlabel("Closeness Centrality")
    plt.ylabel("Final Opinion")
    plt.title("Closeness vs Opinion (r = $(correlations["closeness"]))")
    plt.grid(true)
    
    # Clustering coefficient vs opinion
    plt.subplot(2, 3, 4)
    plt.scatter(clustering_vec, final_opinions, alpha=0.5, color="purple")
    plt.xlabel("Clustering Coefficient")
    plt.ylabel("Final Opinion")
    plt.title("Clustering vs Opinion (r = $(correlations["clustering"])f)")
    plt.grid(true)
    
    # Barra de correlacions
    plt.subplot(2, 3, 5)
    labels = ["Degree", "Betweenness", "Closeness", "Clustering"]
    values = [correlations["degree"], correlations["betweenness"], 
             correlations["closeness"], correlations["clustering"]]
    
    colors = ["blue", "orange", "green", "red"]
    plt.bar(labels, values, color=colors)
    plt.axhline(y=0, linestyle="--", color="gray")
    plt.ylabel("Correlation Coefficient")
    plt.title("Correlation between Structural Properties and Final Opinion")
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig("$PATH_TO_PLOTS/structure_opinion_correlation_$name.png")
    
    # Crear un segon gràfic per visualitzar propietats per comunitat
    plt.figure(figsize=(14, 8))
    
    # Calcular valors mitjans per comunitat
    community_avg_degree = []
    community_avg_clustering = []
    community_avg_opinion = []
    community_sizes = []
    
    for comm in communities
        # Convertir nodes a int per accedir a les opinions
        comm_nodes = [parse(Int, node) for node in comm]
        
        # Calcular la mitjana d'opinions
        avg_opinion = mean(final_opinions[comm_nodes])
        push!(community_avg_opinion, avg_opinion)
        
        # Calcular la mida de la comunitat
        push!(community_sizes, length(comm))
        
        # Calcular la mitjana del grau
        avg_degree = mean([degree_vec[node] for node in comm_nodes])
        push!(community_avg_degree, avg_degree)
        
        # Calcular la mitjana del coeficient d'agrupament
        avg_clustering = mean([clustering_vec[node] for node in comm_nodes])
        push!(community_avg_clustering, avg_clustering)
    end
    
    # Bubble chart: grau vs clustering vs opinió
    plt.scatter(community_avg_degree, community_avg_clustering, 
               s=[size*100 for size in community_sizes], 
               c=community_avg_opinion, cmap="coolwarm", alpha=0.7)
    plt.colorbar(label="Opinió Mitjana")
    plt.xlabel("Community Average Degree")
    plt.ylabel("Community Average Clustering Coefficient")
    plt.title("Community Structure vs Opinion")
    plt.grid(true)
    
    plt.tight_layout()
    plt.savefig("$PATH_TO_PLOTS/community_structure_opinion_$name.png")
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_bridge_regular_dynamics(bridge_opinions, regular_opinions, bridge_variances, regular_variances, name)
    plt.figure(figsize=(12, 10))
    
    # Mean opinion comparison
    plt.subplot(2, 1, 1)
    plt.plot(1:NUM_STEPS, bridge_opinions, label="Nodes Pont", color="red")
    plt.plot(1:NUM_STEPS, regular_opinions, label="Nodes Regulars", color="blue")
    plt.xlabel("Pas Temporal")
    plt.ylabel("Opinió Mitjana")
    plt.title("Evolució de l'Opinió Mitjana: Nodes Pont vs Regulars")
    plt.legend()
    plt.grid(true)
    
    # Variance comparison
    plt.subplot(2, 1, 2)
    plt.plot(1:NUM_STEPS, bridge_variances, label="Bridge Nodes", color="red")
    plt.plot(1:NUM_STEPS, regular_variances, label="Regular Nodes", color="blue")
    plt.xlabel("Step")
    plt.ylabel("Opinion Variance")
    plt.title("Evolution of Opinion Variance: Bridge Nodes vs Regular Nodes")
    plt.legend()
    plt.grid(true)
    
    plt.tight_layout()
    plt.savefig("$PATH_TO_PLOTS/$NETWORK_ANALYSIS_PLOTS/bridge_nodes_analysis_$name.png")
    plt.close()
end
#-----------------------------------------------------------------------------------------------------------------------------------
function visualize_parameter_sensitivity(results)
    # Extreure valors únics dels paràmetres
    λ_values = sort(unique([k[1] for k in keys(results)]))
    trust_values = sort(unique([k[2] for k in keys(results)]))
    glob_period_values = sort(unique([k[3] for k in keys(results)]))
    prob_majority_values = sort(unique([k[4] for k in keys(results)]))
    
    # Crear directori si no existeix
    plots_dir = "$PATH_TO_PLOTS/parameter_sensitivity"
    mkpath(plots_dir)
    
    # -------------------------------------------------------------------------
    # 1. Heatmap: Efecte de λ i TRUST_LEVEL en la polarització
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 8))
    
    # Crear matriu per al heatmap amb valors per defecte dels altres paràmetres
    default_period = 500  # Període global típic
    default_prob = 0.7    # Conformitat inicial típica
    
    # Trobar els índexs més propers als valors per defecte
    period_idx = findmin(abs.(glob_period_values .- default_period))[2]
    prob_idx = findmin(abs.(prob_majority_values .- default_prob))[2]
    
    global default_period = glob_period_values[period_idx]
    global default_prob = prob_majority_values[prob_idx]
    
    # Crear la matriu de dades
    polarization_matrix = zeros(length(λ_values), length(trust_values))
    
    for (i, λ_val) in enumerate(λ_values)
        for (j, trust_val) in enumerate(trust_values)
            key = (λ_val, trust_val, default_period, default_prob)
            if haskey(results, key)
                polarization_matrix[i, j] = results[key]["polarization"]
            end
        end
    end
    
    # Crear heatmap
    im = plt.imshow(polarization_matrix, cmap="viridis", aspect="auto", origin="lower")
    cbar = plt.colorbar(im, label="Índex de Polarització")
    
    plt.title("Efecte de λ i TRUST_LEVEL en la Polarització\n" *
              "(Període Global: $default_period, Conformitat Inicial: $(round(default_prob*100))%)")
    plt.xlabel("TRUST_LEVEL (llindar de confiança)")
    plt.ylabel("λ (sensibilitat a distància d'opinió)")
    
    # Afegir anotacions explicatives
    plt.figtext(0.02, 0.02, "λ baix: Menys sensible a diferències d'opinió\nλ alt: Més sensible a diferències d'opinió", 
                bbox=Dict("facecolor" => "white", "alpha" => 0.7), fontsize=9)
    plt.figtext(0.65, 0.02, "TRUST_LEVEL baix: Més fàcil canviar d'opinió\nTRUST_LEVEL alt: Més difícil canviar d'opinió", 
                bbox=Dict("facecolor" => "white", "alpha" => 0.7), fontsize=9)
    
    plt.tight_layout()
    plt.savefig("$plots_dir/lambda_trust_polarization.png", dpi=300)
    plt.close()
    
    # -------------------------------------------------------------------------
    # 2. Efecte de la Conformitat Inicial en la Polarització
    # -------------------------------------------------------------------------
    plt.figure(figsize=(12, 7))
    
    # Seleccionar un subconjunt representatiu de valors de λ si n'hi ha molts
    λ_subset = if length(λ_values) <= 5
        λ_values
    else
        # Seleccionar aproximadament 5 valors distribuïts uniformement
        step = max(1, div(length(λ_values), 5))
        λ_values[1:step:end]
    end
    
    for λ_val in λ_subset
        polarization_values = []
        
        for prob in prob_majority_values
            # Valors per defecte per als altres paràmetres
            default_trust = 0.7  # Valor comú de TRUST_LEVEL
            
            # Trobar índex més proper si no existeix el valor exacte
            trust_idx = findmin(abs.(trust_values .- default_trust))[2]
            global default_trust = trust_values[trust_idx]
            
            key = (λ_val, default_trust, default_period, prob)
            if haskey(results, key)
                push!(polarization_values, results[key]["polarization"])
            else
                push!(polarization_values, NaN)
            end
        end
        
        plt.plot(prob_majority_values, polarization_values, marker="o", 
                linewidth=2, label="λ=$(round(λ_val, digits=2))")
    end
    
    plt.title("Efecte de la Conformitat Inicial en la Polarització")
    plt.xlabel("Conformitat Inicial (% nodes amb opinió majoritària)")
    plt.ylabel("Índex de Polarització")
    plt.grid(true, linestyle="--", alpha=0.7)
    plt.legend(title="Sensibilitat (λ)")
    
    # Ajustar etiquetes de l'eix X per mostrar percentatges
    plt.xticks(prob_majority_values, ["$(round(p*100))%" for p in prob_majority_values])
    
    # Afegir anotació explicativa
    plt.figtext(0.02, 0.02, "Conformitat alta: Més nodes inicien amb la mateixa opinió\n" *
                          "Conformitat baixa: Més diversitat d'opinions inicials", 
                bbox=Dict("facecolor" => "white", "alpha" => 0.7), fontsize=9)
    
    plt.tight_layout()
    plt.savefig("$plots_dir/conformity_polarization.png", dpi=300)
    plt.close()
    
    # -------------------------------------------------------------------------
    # 3. Efecte del Període Global en el Temps de Convergència
    # -------------------------------------------------------------------------
    plt.figure(figsize=(12, 7))
    
    for λ_val in λ_subset
        convergence_times = []
        
        for period in glob_period_values
            # Valors per defecte
            key = (λ_val, default_trust, period, default_prob)
            if haskey(results, key)
                push!(convergence_times, results[key]["convergence_step"])
            else
                push!(convergence_times, NaN)
            end
        end
        
        plt.plot(glob_period_values, convergence_times, marker="s", 
                linewidth=2, label="λ=$(round(λ_val, digits=2))")
    end
    
    plt.title("Efecte del Període d'Influència Global en el Temps de Convergència")
    plt.xlabel("Període d'Influència Global (passos)")
    plt.ylabel("Passos fins a Convergència")
    plt.grid(true, linestyle="--", alpha=0.7)
    plt.legend(title="Sensibilitat (λ)")
    
    # Afegir línia vertical en període = 0 (sense influència global)
    if 0 in glob_period_values
        plt.axvline(x=0, color="red", linestyle="--", alpha=0.5)
        plt.text(5, plt.ylim()[1] + 0.05*(plt.ylim()[1]-plt.ylim()[0]), 
                "Sense influència global", color="red", rotation=90)
    end
    
    # Afegir anotació explicativa
    plt.figtext(0.02, 0.02, "Període baix: Influència global freqüent\n" *
                          "Període alt: Influència global poc freqüent\n" *
                          "Període 0: Sense influència global", 
                bbox=Dict("facecolor" => "white", "alpha" => 0.7), fontsize=9)
    
    plt.tight_layout()
    plt.savefig("$plots_dir/global_period_convergence.png", dpi=300)
    plt.close()
    
    # -------------------------------------------------------------------------
    # 4. Resum de Paràmetres Crítics
    # -------------------------------------------------------------------------
    # Trobar valors de paràmetres que produeixen la polarització màxima i mínima
    max_polarization = -Inf
    min_polarization = Inf
    max_key = nothing
    min_key = nothing
    
    for key in keys(results)
        if haskey(results[key], "polarization")
            pol = results[key]["polarization"]
            if pol > max_polarization
                max_polarization = pol
                max_key = key
            end
            if pol < min_polarization
                min_polarization = pol
                min_key = key
            end
        end
    end
    
    # Crear una taula visual amb els paràmetres crítics
    plt.figure(figsize=(10, 6))
    plt.axis("off")
    
    if max_key !== nothing && min_key !== nothing
        title_text = "Paràmetres Crítics per a la Polarització"
        max_text = "Màxima Polarització ($(round(max_polarization, digits=3))):\n" *
                   "  λ = $(max_key[1])\n" *
                   "  TRUST_LEVEL = $(max_key[2])\n" *
                   "  Període Global = $(max_key[3])\n" *
                   "  Conformitat Inicial = $(round(max_key[4]*100))%"
        
        min_text = "Mínima Polarització ($(round(min_polarization, digits=3))):\n" *
                   "  λ = $(min_key[1])\n" *
                   "  TRUST_LEVEL = $(min_key[2])\n" *
                   "  Període Global = $(min_key[3])\n" *
                   "  Conformitat Inicial = $(round(min_key[4]*100))%"
        
        params_text = "Sumari dels Paràmetres:\n\n" *
                      "λ: Controla la sensibilitat a la distància d'opinió\n" *
                      "   P_i = 1/(1 + e^((1 - d_oi)/λ))\n\n" *
                      "TRUST_LEVEL: Llindar mínim per adoptar opinions\n" *
                      "   Si Pi > TRUST_LEVEL → el node adopta l'opinió\n\n" *
                      "Període Global: Freqüència d'influència global\n" *
                      "   Cada X passos s'aplica una influència externa\n\n" *
                      "Conformitat Inicial: Percentatge de nodes que comencen\n" *
                      "   amb l'opinió majoritària de la seva comunitat"
        
        plt.text(0.5, 0.95, title_text, fontsize=16, fontweight="bold", ha="center")
        plt.text(0.05, 0.75, max_text, fontsize=12, bbox=Dict("facecolor" => "lightyellow", "alpha" => 0.5))
        plt.text(0.05, 0.55, min_text, fontsize=12, bbox=Dict("facecolor" => "lightblue", "alpha" => 0.5))
        plt.text(0.05, 0.2, params_text, fontsize=10, va="top")
    else
        plt.text(0.5, 0.5, "No hi ha dades suficients per trobar paràmetres crítics", ha="center")
    end
    
    plt.tight_layout()
    plt.savefig("$plots_dir/critical_parameters.png", dpi=300)
    plt.close()
    
    println("Anàlisi de sensibilitat de paràmetres completada. Gràfics desats a: $plots_dir")
    
    return nothing
end
#-----------------------------------------------------------------------------------------------------------------------------------
end