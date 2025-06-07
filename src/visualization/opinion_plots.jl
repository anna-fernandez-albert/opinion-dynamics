module Visualization

include("../constants.jl")
using PyCall
using Plots
using Dierckx
using Statistics
using LinearAlgebra
using Polynomials
@pyimport networkx as nx
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.cm as cm
@pyimport pyvis.network as pyvis

export plot_network, plot_network_stable_nodes, plot_opinion_evolution, plot_correlation_bars, heatmap_consensus_lambda_tau, plot_execution_vs_T, plot_fluctuant_vs_tau
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
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_network_instable_state(graph, analysis, path)
    println("Plotting network with instable nodes...")
    G = pyvis.Network(height="1000px", width="100%")
    pos = nx.spring_layout(graph)
    colors = Dict(-1 => "red", 0 => "blue", 1 => "green")
    for (node, opinion) in analysis
        G.add_node(node, x=pos[node][1]*10000, y=pos[node][2]*10000, size=150, physics=false, color=colors[opinion])
    end
    for (u, v) in graph.edges
        G.add_edge(u, v, physics=false, width=0.5, color="gray")     
    end
    println("Saving network visualization to: $path")
    G.save_graph(path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_opinion_evolution(opinion_history, path)
    avg_opinion = abs.(mean(opinion_history, dims=2))
    p = plot(vec(avg_opinion), xlabel="Time Step", ylabel="Mean Opinion")
    savefig(p, path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
function plot_correlation_bars(correlations, path)
    parameter_names = collect(keys(correlations))
    values = [correlations[k] for k in parameter_names]
    p = bar(parameter_names, values, xlabel="Parameter", ylabel="Correlation")
    savefig(p, path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
function heatmap_consensus_lambda_tau(results, path)
    # if the key trust_value does not exist, change trust_level to trust_value
    if !haskey(results[1], "trust_value")
        for r in results
            r["trust_value"] = r["trust_level"]
        end
    end
    λ_vals = unique([r["λ"] for r in results])
    trust_vals = unique([r["trust_value"] for r in results])
    heatmap_data = [mean([r["consensus"] for r in results if r["λ"] == λ && r["trust_value"] == τ]) for τ in trust_vals, λ in λ_vals]
    p = heatmap(λ_vals, trust_vals, heatmap_data, xlabel="λ", ylabel="τ", size=(620, 600), color=:plasma, aspect_ratio=:auto)
    savefig(p, path)
end
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_execution_vs_T(results, path)
    # Extract all data points
    T_vals = [r["global_period_value"] for r in results]
    t_exec_vals = [r["t_execution"] for r in results]
    
    p = plot(xlabel="Global Period T", ylabel="Execution Time (steps)")
    # mean trend line
    T_unique = sort(unique(T_vals))
    mean_times = [mean([t_exec_vals[i] for i in 1:length(T_vals) if T_vals[i] == T]) for T in T_unique]
    
    plot!(T_unique, mean_times, color=:darkblue, linewidth=1.5, label="Mean Trend", lineestyle=:dash)
    
    savefig(p, path)
end
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_fluctuant_vs_tau(results, path)
    # Plot the mean fluctuation fraction against the trust value
    grouped = Dict{Float64, Vector{Float64}}()
    for r in results
        τ = r["trust_value"]
        push!(get!(grouped, τ, []), r["flip_fraction"])
    end

    # Extract data
    taus = sort(collect(keys(grouped)))
    mean_fluctuations = [mean(grouped[τ]) for τ in taus]
    std_fluctuations = [std(grouped[τ]) for τ in taus]

    # Add confidence interval
    #p = plot(taus, mean_fluctuations .+ std_fluctuations, fillto=mean_fluctuations .- std_fluctuations, alpha=0.3, color=:steelblue, label="±σ")
    
    p = plot(taus, mean_fluctuations, xlabel="Trust Value τ", ylabel="Mean Fluctuant Fraction", color=:steelblue, linewidth=2)
    plot!(p, legend=false)
    savefig(p, path)
end
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_fluctuant_vs_lambda(results, path)
    # Group by lambda value
    grouped = Dict{Float64, Vector{Float64}}()
    for r in results
        λ = r["λ"]
        push!(get!(grouped, λ, []), r["flip_fraction"])
    end
    
    # Extract data
    λs = sort(collect(keys(grouped)))
    mean_fluctuations = [mean(grouped[λ]) for λ in λs]
    std_fluctuations = [std(grouped[λ]) for λ in λs]
    # Add confidence interval
    #p = plot(λs, mean_fluctuations .+ std_fluctuations, fillto=mean_fluctuations .- std_fluctuations, alpha=0.3, color=:crimson, label="±σ")
    
    p = plot(λs, mean_fluctuations, xlabel="λ", ylabel="Mean Fluctuant Fraction", color=:crimson, linewidth=2)

    plot!(p, legend=false)
    xlims!(p, (0.2, maximum(λs)))
    savefig(p, path)
end
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_consensus_vs_tau(results, path)
    # Group by trust value
    grouped = Dict{Float64, Vector{Float64}}()
    for r in results
        τ = r["trust_value"]
        push!(get!(grouped, τ, []), r["consensus"])
    end

    # Extract data
    taus = sort(collect(keys(grouped)))
    mean_consensus = [mean(grouped[τ]) for τ in taus]
    std_consensus = [std(grouped[τ]) for τ in taus]

    # Mean line with confidence interval
    p = plot(taus, mean_consensus .+ std_consensus, fillto=mean_consensus .- std_consensus, alpha=0.3, color=:steelblue, label="±σ")
    
    # Mean trend line
    plot!(p, taus, mean_consensus, color=:steelblue, linewidth=2, label="Mean Consensus")

    plot!(p, xlabel="Trust Value τ", ylabel="Mean Consensus")
    plot!(p, legend=:bottomright)
    ylims!(p, (0.25, 1.25))

    savefig(p, path)
end
#--Only for Parameter Sensibility---------------------------------------------------------------------------------------------------
function plot_consensus_vs_lambda(results, path)
    # Group by lambda value
    grouped = Dict{Float64, Vector{Float64}}()
    for r in results
        λ = r["λ"]
        push!(get!(grouped, λ, []), r["consensus"])
    end
    
    # Extract data
    λs = sort(collect(keys(grouped)))
    mean_consensus = [mean(grouped[λ]) for λ in λs]
    std_consensus = [std(grouped[λ]) for λ in λs]
    
    # Confidence interval
    p = plot(λs, mean_consensus .+ std_consensus, fillto=mean_consensus .- std_consensus, alpha=0.3, color=:crimson, label="±σ")
    
    # Mean trend line (the key result!)
    plot!(p, λs, mean_consensus, color=:crimson, linewidth=2, label="Mean Consensus")

    # Add critical point annotation
    # Compute the second derivative to find the inflection point
    λ_critical = λs[argmax(diff(mean_consensus))]
    vline!(p, [λ_critical], color=:black, linestyle=:dash, linewidth=1, alpha=0.8, label="Critical Point (λ = $λ_critical)")
    
    # Formatting
    plot!(p, legend=:bottomright, grid=true, gridwidth=1, gridcolor=:lightgray)
    ylims!(p, (0.25, 1.25))
    xlims!(p, (minimum(λs), maximum(λs)))
    xlabel!(p, "λ")
    ylabel!(p, "Mean Consensus")
    savefig(p, path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_assortativity(results, path)
    """
    Plot Assortativity vs Consensus.
    Each point is a network (N, μ, k_avg) with its assortativity and mean consensus (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its assortativity and mean consensus
    assortativities = Float64[]
    consensus_means = Float64[]
    k_avg_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(assortativities, first(data)["assortativity"])
        push!(consensus_means, mean([r["consensus"] for r in data]))
        push!(k_avg_vals, k_avg)
    end

    # Plot
    p = scatter(assortativities, consensus_means, zcolor=k_avg_vals,
                xlabel="Assortativity (r)", ylabel="Mean Consensus",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="k average", label="")

    # Trend line
    trend_x = collect(minimum(assortativities):0.01:maximum(assortativities))
    pf = fit(assortativities, consensus_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_std_dev_community_size(results, path)
    """
    Plot Standard Deviation of Community Size vs Consensus.
    Each point is a network (N, μ, k_avg) with its standard deviation of community size and mean consensus (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its std_dev_community_size and mean consensus
    std_dev_sizes = Float64[]
    consensus_means = Float64[]
    mu_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(std_dev_sizes, first(data)["std_dev_community_size"])
        push!(consensus_means, mean([r["consensus"] for r in data]))
        push!(mu_vals, μ)
    end

    # Plot
    p = scatter(std_dev_sizes, consensus_means, zcolor=mu_vals,
                xlabel="Standard Deviation of Community Size", ylabel="Mean Consensus",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="μ (Mixing)", label="")

    # Trend line
    trend_x = collect(minimum(std_dev_sizes):0.01:maximum(std_dev_sizes))
    pf = fit(std_dev_sizes, consensus_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_clustering_coefficient(results, path)
    """
    Plot Clustering Coefficient vs Consensus.
    Each point is a network (N, μ, k_avg) with its clustering coefficient and mean consensus (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its clustering coefficient and mean consensus
    clustering_coeffs = Float64[]
    consensus_means = Float64[]
    mu_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(clustering_coeffs, first(data)["clustering_coefficient"])
        push!(consensus_means, mean([r["consensus"] for r in data]))
        push!(mu_vals, μ)
    end

    # Plot
    p = scatter(clustering_coeffs, consensus_means, zcolor=mu_vals,
                xlabel="Clustering Coefficient", ylabel="Mean Consensus",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="μ (Mixing)", label="")

    # Trend line
    trend_x = collect(minimum(clustering_coeffs):0.01:maximum(clustering_coeffs))
    pf = fit(clustering_coeffs, consensus_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_density_consensus_heatmap(results, path)
    Ns = sort(unique([r["N"] for r in results]))
    k_vals = sort(unique([r["k_avg"] for r in results]))
    
    # Create heatmap data with proper handling of missing combinations
    heatmap_data = zeros(length(k_vals), length(Ns))
    for (i, k) in enumerate(k_vals)
        for (j, N) in enumerate(Ns)
            matching_results = filter(r -> r["N"] == N && r["k_avg"] == k, results)
            if !isempty(matching_results)
                heatmap_data[i, j] = mean([r["consensus"] for r in matching_results])
            else
                heatmap_data[i, j] = NaN  # Handle missing data
            end
        end
    end
    
    p = heatmap(Ns, k_vals, heatmap_data, xlabel="Network Size (N)",  ylabel="Average Degree ⟨k⟩", size=(625, 600), color=:plasma, aspect_ratio=:auto)
    
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_polarization_vs_modularity(results, path)
    """
    Plot Polarization Between Communities vs Modularity.
    Each point is a network (N, μ, k_avg) with its modularity and mean polarization (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its modularity and mean polarization
    modularities = Float64[]
    polarization_means = Float64[]
    mu_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(modularities, first(data)["modularity"])
        push!(polarization_means, mean([r["polarization_between_communities"] for r in data]))
        push!(mu_vals, μ)
    end

    # Plot
    p = scatter(modularities, polarization_means, zcolor=mu_vals,
                xlabel="Modularity", ylabel="Mean Polarization Between Communities",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="μ (Mixing)")

    # Trend line
    trend_x = collect(minimum(modularities):0.01:maximum(modularities))
    pf = fit(modularities, polarization_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:bottomright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_modularity(results, path)
    """
    Plot Consensus vs Modularity.
    Each point is a network (N, μ, k_avg) with its modularity and mean consensus (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its modularity and mean consensus
    modularities = Float64[]
    consensus_means = Float64[]
    mu_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(modularities, first(data)["modularity"])
        push!(consensus_means, mean([r["consensus"] for r in data]))
        push!(mu_vals, μ)
    end

    # Plot
    p = scatter(modularities, consensus_means, zcolor=mu_vals,
                xlabel="Modularity", ylabel="Mean Consensus",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="μ (Mixing)", label="")

    # Trend line
    trend_x = collect(minimum(modularities):0.01:maximum(modularities))
    pf = fit(modularities, consensus_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_av_path_length(results, path)
    """
    Plot Consensus vs Average Path Length.
    Each point is a network (N, μ, k_avg) with its average path length and mean consensus (averaged over λ and τ).
    """
    # Group results by network configuration (N, μ, k_avg)
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # For each network, get its average path length and mean consensus
    av_path_lengths = Float64[]
    consensus_means = Float64[]
    k_avg_vals = Float64[]
    for ((N, μ, k_avg), data) in grouped
        push!(av_path_lengths, first(data)["av_path_length"])
        push!(consensus_means, mean([r["consensus"] for r in data]))
        push!(k_avg_vals, k_avg)
    end

    # Plot
    p = scatter(av_path_lengths, consensus_means, zcolor=k_avg_vals,
                xlabel="Average Path Length", ylabel="Mean Consensus",
                markersize=6, alpha=0.8, color=:viridis, colorbar_title="k average", label="")

    # Trend line
    trend_x = collect(minimum(av_path_lengths):0.01:maximum(av_path_lengths))
    pf = fit(av_path_lengths, consensus_means, 1)
    trend_y = pf.(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    plot!(p, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_trust_level_different_net(results, path)
    """Trust effect across different network structures (mean consensus over all λ for each τ)"""
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # Generate a color gradient from light red to dark red
    n_lines = length(grouped)
    reds = cgrad(:reds, n_lines, rev=true)  # from light to dark red

    p = plot(xlabel="Trust Level τ", ylabel="Consensus")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        # Group by trust_level and average consensus over all λ
        trust_groups = Dict{Float64, Vector{Float64}}()
        for r in data
            τ = r["trust_level"]
            push!(get!(trust_groups, τ, []), r["consensus"])
        end
        trust_levels = sort(collect(keys(trust_groups)))
        consensus_means = [mean(trust_groups[τ]) for τ in trust_levels]
        plot!(p, trust_levels, consensus_means, label="N=$N, μ=$μ, k_avg=$k_avg", linewidth=1.5, alpha=0.8, color=reds[i])
        i += 1
    end
    plot!(p, legend=false)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_lambda_different_net(results, path)
    """Lambda phase transitions across network structures (mean consensus over all τ for each λ)"""
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # Generate a color gradient from light blue to dark blue
    n_lines = length(grouped)
    blues = cgrad(:blues, n_lines, rev=true)  # from light to dark blue

    p = plot(xlabel="λ", ylabel="Consensus")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        # Group by λ and average consensus over all τ
        lambda_groups = Dict{Float64, Vector{Float64}}()
        for r in data
            λ = r["λ"]
            push!(get!(lambda_groups, λ, []), r["consensus"])
        end
        lambdas = sort(collect(keys(lambda_groups)))
        consensus_means = [mean(lambda_groups[λ]) for λ in lambdas]
        plot!(p, lambdas, consensus_means, label="N=$N, μ=$μ, k_avg=$k_avg", linewidth=1.5, alpha=0.8, color=blues[i])
        i += 1
    end
    plot!(p, legend=false)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_fluctuant_vs_lambda_different_net(results, path)
    """Instability peaks across network structures (mean fluctuant fraction over all τ for each λ)"""
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # Generate a color gradient from light blue to dark blue
    n_lines = length(grouped)
    blues = cgrad(:blues, n_lines, rev=true)  # from light to dark blue

    p = plot(xlabel="λ", ylabel="Fluctuant Fraction")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        # Group by λ and average flip_fraction over all τ
        lambda_groups = Dict{Float64, Vector{Float64}}()
        for r in data
            λ = r["λ"]
            push!(get!(lambda_groups, λ, []), r["flip_fraction"])
        end
        lambdas = sort(collect(keys(lambda_groups)))
        fluctuant_means = [mean(lambda_groups[λ]) for λ in lambdas]
        plot!(p, lambdas, fluctuant_means, label="N=$N, μ=$μ, k_avg=$k_avg", linewidth=1.5, alpha=0.8, color=blues[i])
        i += 1
    end
    plot!(p, legend=false)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------    
function plot_fluctuant_vs_trust_level_different_net(results, path)
    """Trust effect on system stability across structures (mean fluctuant fraction over all λ for each τ)"""
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Dict}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), r)
    end

    # Generate a color gradient from light red to dark red
    n_lines = length(grouped)
    reds = cgrad(:reds, n_lines, rev=true)  # from light to dark red

    p = plot(xlabel="Trust Level τ", ylabel="Fluctuant Fraction")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        # Group by trust_level and average flip_fraction over all λ
        trust_groups = Dict{Float64, Vector{Float64}}()
        for r in data
            τ = r["trust_level"]
            push!(get!(trust_groups, τ, []), r["flip_fraction"])
        end
        trust_levels = sort(collect(keys(trust_groups)))
        fluctuant_means = [mean(trust_groups[τ]) for τ in trust_levels]
        plot!(p, trust_levels, fluctuant_means, label="N=$N, μ=$μ, k_avg=$k_avg", linewidth=1.5, alpha=0.8, color=reds[i])
        i += 1
    end
    plot!(p, legend=false)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_lambda_modularity(results, path)
    """
    Lambda phase transitions across modularity types (mean consensus over all τ for each λ).
    Only three lines are plotted: one for each modularity type (low, medium, high), each being the mean over all networks of that type.
    Modularity types are determined by the value of 'modularity' in the results:
        low modularity: modularity < 0.2
        medium modularity: 0.2 <= modularity < 0.4
        high modularity: modularity >= 0.4
    """
    # Group results by modularity type and λ
    consensus_dict = Dict("low" => Dict{Float64, Vector{Float64}}(),
                         "medium" => Dict{Float64, Vector{Float64}}(),
                         "high" => Dict{Float64, Vector{Float64}}())

    for r in results
        m = r["modularity"]
        λ = r["λ"]
        consensus = r["consensus"]
        if m < 0.2
            t = "low"
        elseif m < 0.4
            t = "medium"
        else
            t = "high"
        end
        push!(get!(consensus_dict[t], λ, Float64[]), consensus)
    end

    # For each modularity type, compute mean consensus for each λ
    plot_colors = [:lightblue, :dodgerblue, :navy]
    labels = ["Low Modularity", "Medium Modularity", "High Modularity"]
    p = plot(xlabel="λ", ylabel="Consensus")
    for (i, t) in enumerate(["low", "medium", "high"])
        λs = sort(collect(keys(consensus_dict[t])))
        means = [mean(consensus_dict[t][λ]) for λ in λs]
        # Add critical point annotation
        λ_critical = λs[argmax(diff(means))]
        vline!(p, [λ_critical], color=:black, linestyle=:dash, linewidth=1, alpha=0.8, label="Critical Point (λ = $λ_critical)")
        plot!(p, λs, means, label=labels[i], linewidth=1.5, alpha=0.8, color=plot_colors[i])
    end
    plot!(p, legend=:bottomright)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_trust_level_modularity(results, path)
    """
    Trust effect on consensus across modularity types (mean consensus over all λ for each τ).
    Only three lines are plotted: one for each modularity type (low, medium, high), each being the mean over all networks of that type.
    Modularity types are determined by the value of 'modularity' in the results:
        low modularity: modularity < 0.2
        medium modularity: 0.2 <= modularity < 0.4
        high modularity: modularity >= 0.4
    """
    # Group results by modularity type and trust_value
    consensus_dict = Dict("low" => Dict{Float64, Vector{Float64}}(),
                         "medium" => Dict{Float64, Vector{Float64}}(),
                         "high" => Dict{Float64, Vector{Float64}}())

    for r in results
        m = r["modularity"]
        τ = r["trust_value"]
        consensus = r["consensus"]
        if m < 0.2
            t = "low"
        elseif m < 0.4
            t = "medium"
        else
            t = "high"
        end
        push!(get!(consensus_dict[t], τ, Float64[]), consensus)
    end

    # For each modularity type, compute mean consensus for each τ
    plot_colors = [:lightcoral, :indianred, :firebrick]
    labels = ["Low Modularity", "Medium Modularity", "High Modularity"]
    p = plot(xlabel="Trust Level τ", ylabel="Consensus")
    for (i, t) in enumerate(["low", "medium", "high"])
        τs = sort(collect(keys(consensus_dict[t])))
        means = [mean(consensus_dict[t][τ]) for τ in τs]
        # Add the line tendency and its function to the plot
        plot!(p, τs, means, label=labels[i], linewidth=1.5, alpha=0.8, color=plot_colors[i])
        # Add a trend line
        trend_x = collect(minimum(τs):0.01:maximum(τs))
        pf = fit(τs, means, 1)
        trend_y = pf.(trend_x)
        plot!(p, trend_x, trend_y, color=:black, linestyle=:dash, linewidth=1.5, alpha=0.8, label="")
        # Add the function of the trend line
        trend_func = pf
        coeffs_vec = coeffs(trend_func)
        trend_eq = "y = $(round(coeffs_vec[2], digits=3)) * x + $(round(coeffs_vec[1], digits=3))"
        annotate!(p, [(mean(τs), mean(means), text(trend_eq, :left, 10, :black))])
    end
    plot!(p, legend=:bottomright)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------    
function plot_fluctuant_vs_lambda_modularity(results, path)
    """
    Instability peaks across modularity types (mean fluctuant fraction over all τ for each λ).
    Only three lines are plotted: one for each modularity type (low, medium, high), each being the mean over all networks of that type.
    Modularity types are determined by the value of 'modularity' in the results:
        low modularity: modularity < 0.2
        medium modularity: 0.2 <= modularity < 0.4
        high modularity: modularity >= 0.4
    """
    # Group results by modularity type and λ
    fluctuant_dict = Dict("low" => Dict{Float64, Vector{Float64}}(),
                         "medium" => Dict{Float64, Vector{Float64}}(),
                         "high" => Dict{Float64, Vector{Float64}}())

    for r in results
        m = r["modularity"]
        λ = r["λ"]
        flip = r["flip_fraction"]
        if m < 0.2
            t = "low"
        elseif m < 0.4
            t = "medium"
        else
            t = "high"
        end
        push!(get!(fluctuant_dict[t], λ, Float64[]), flip)
    end

    # For each modularity type, compute mean fluctuant fraction for each λ
    plot_colors = [:lightblue, :dodgerblue, :navy]
    labels = ["Low Modularity", "Medium Modularity", "High Modularity"]
    p = plot(xlabel="λ", ylabel="Fluctuant Fraction")
    for (i, t) in enumerate(["low", "medium", "high"])
        λs = sort(collect(keys(fluctuant_dict[t])))
        means = [mean(fluctuant_dict[t][λ]) for λ in λs]
        plot!(p, λs, means, label=labels[i], linewidth=1.5, alpha=0.8, color=plot_colors[i])
    end
    plot!(p, legend=:topleft)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_fluctuant_vs_trust_level_modularity(results, path)
    """
    Trust effect on system stability across modularity types (mean fluctuant fraction over all λ for each τ).
    Only three lines are plotted: one for each modularity type (low, medium, high), each being the mean over all networks of that type.
    Modularity types are determined by the value of 'modularity' in the results:
        low modularity: modularity < 0.2
        medium modularity: 0.2 <= modularity < 0.4
        high modularity: modularity >= 0.4
    """
    # Group results by modularity type and trust_value
    fluctuant_dict = Dict("low" => Dict{Float64, Vector{Float64}}(),
                         "medium" => Dict{Float64, Vector{Float64}}(),
                         "high" => Dict{Float64, Vector{Float64}}())

    for r in results
        m = r["modularity"]
        τ = r["trust_value"]
        flip = r["flip_fraction"]
        if m < 0.2
            t = "low"
        elseif m < 0.4
            t = "medium"
        else
            t = "high"
        end
        push!(get!(fluctuant_dict[t], τ, Float64[]), flip)
    end

    # For each modularity type, compute mean fluctuant fraction for each τ
    plot_colors = [:lightcoral, :indianred, :firebrick]
    labels = ["Low Modularity", "Medium Modularity", "High Modularity"]
    p = plot(xlabel="Trust Level τ", ylabel="Fluctuant Fraction")
    for (i, t) in enumerate(["low", "medium", "high"])
        τs = sort(collect(keys(fluctuant_dict[t])))
        means = [mean(fluctuant_dict[t][τ]) for τ in τs]
        plot!(p, τs, means, label=labels[i], linewidth=1.5, alpha=0.8, color=plot_colors[i])
    end
    plot!(p, legend=:topleft)
    savefig(p, path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
end