module Visualization

include("../constants.jl")
using PyCall
using Plots
using Dierckx
using Statistics
using LinearAlgebra
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
function plot_consensus_vs_mu(results, path)
    """Consensus vs Mixing Parameter μ for different network configurations"""
    
    # Group by network configuration
    grouped = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}()
    for r in results
        k = (r["N"], r["k_avg"])
        push!(get!(grouped, k, []), (r["μ"], r["consensus"]))
    end

    colors = [:steelblue, :crimson, :forestgreen, :darkorange, :purple, :brown]
    
    p = plot(xlabel="μ (Mixing Parameter)", ylabel="Consensus")
    color_idx = 1
    for ((N, k_avg), data) in grouped
        sorted_data = sort(data, by=x->x[1])
        mus = [x[1] for x in sorted_data]
        consensus_vals = [x[2] for x in sorted_data]
        
        plot!(p, mus, consensus_vals, label="", color=colors[color_idx % length(colors) + 1], linewidth=1, alpha=0.8)
        scatter!(p, mus, consensus_vals, color=colors[color_idx % length(colors) + 1], markersize=3, alpha=0.7, label="")
        
        color_idx += 1
    end
    
    plot!(p, legend=:false)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_assortativity(results, path)
    """Consensus vs Assortativity with trend analysis"""
    
    # Sort and extract data
    results = sort(results, by=r -> r["assortativity"])
    assortativities = [r["assortativity"] for r in results]
    consensus_vals = [r["consensus"] for r in results]
    
    # Color points by μ for additional insight
    mu_vals = [r["μ"] for r in results]
    
    p = scatter(assortativities, consensus_vals, zcolor=mu_vals,
               xlabel="Assortativity (r)", ylabel="Consensus",
               markersize=4, alpha=0.7, color=:viridis, colorbar_title="μ (Mixing)")
    
    # Add trend line
    trend_x = collect(minimum(assortativities):0.01:maximum(assortativities))
    trend_y = polyfit(assortativities, consensus_vals, 1)(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)

    # Add critical point annotation
    critical_idx = argmax(diff(trend_y))
    print("Critical point (consensus vs assortativity) at r = $(trend_x[critical_idx])")
    
    plot!(p, legend=:bottomright)
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
    results = sort(results, by=r -> r["modularity"])
    modularities = [r["modularity"] for r in results]
    polar_vals = [r["polarization_between_communities"] for r in results]

    p = scatter(modularities, polar_vals, xlabel="Modularity", ylabel="Polarization Between Communities",
               markersize=2, alpha=0.7)

    # Add mean trend line
    trend_x = collect(minimum(modularities):0.01:maximum(modularities))
    trend_y = polyfit(modularities, polar_vals, 1)(trend_x)
    plot!(p, trend_x, trend_y, color=:red, linewidth=1, label="Trend Line", linestyle=:dash)

    p = scatter(modularities, polar_vals, xlabel="Modularity", ylabel="Polarization between Communities")
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_av_path_length(results, path)
    """Consensus vs Average Path Length with network size information"""
    
    results = sort(results, by=r -> r["av_path_length"])
    av_path_lengths = [r["av_path_length"] for r in results]
    consensus_vals = [r["consensus"] for r in results]
    
    p = scatter(av_path_lengths, consensus_vals, xlabel="Average Path Length", ylabel="Consensus", markersize=4, alpha=0.7, color=:viridis)
    
    # Add trend line
    trend_x = collect(minimum(av_path_lengths):0.01:maximum(av_path_lengths))
    trend_y = polyfit(av_path_lengths, consensus_vals, 1)(trend_x)
    
    plot!(p, trend_x, trend_y, color=:red, linewidth=2, label="Trend Line", linestyle=:dash)
    
    plot!(p, grid=true, gridwidth=1, gridcolor=:lightgray, legend=:topright)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_trust_level_different_net(results, path)
    """Trust effect across different network structures"""
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Tuple{Float64, Float64}}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), (r["trust_level"], r["consensus"]))
    end

    # Generate a color gradient from light red to dark red
    n_lines = length(grouped)
    reds = cgrad(:reds, n_lines, rev=true)  # from light to dark red

    p = plot(xlabel="Trust Level τ", ylabel="Consensus")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        sorted_data = sort(data, by=x->x[1])  # sort by trust level
        trust_levels = [x[1] for x in sorted_data]
        consensus_vals = [x[2] for x in sorted_data]
        plot!(p, trust_levels, consensus_vals, label="N=$N, μ=$μ, k_avg=$k_avg", linewidth=1.5, alpha=0.8, color=reds[i])
        i += 1
    end
    plot!(p, legend=false)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_consensus_vs_lambda_different_net(results, path)
    """Lambda phase transitions across network structures"""
    
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Tuple{Float64, Float64}}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), (r["λ"], r["consensus"]))
    end

    # Generate a color gradient from light blue to dark blue
    n_lines = length(grouped)
    blues = cgrad(:blues, n_lines, rev=true)  # from light to dark blue

    p = plot(xlabel="λ", ylabel="Consensus")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        sorted_data = sort(data, by=x->x[1])  # sort by λ
        lambdas = [x[1] for x in sorted_data]
        consensus_vals = [x[2] for x in sorted_data]
        plot!(p, lambdas, consensus_vals, label="N=$N, μ=$μ, k_avg=$k_avg", color=blues[i], linewidth=1.5, alpha=0.8)
        i += 1
    end
    plot!(p, legend=false)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------
function plot_fluctuant_vs_lambda_different_net(results, path)
    """Instability peaks across network structures"""
    
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Tuple{Float64, Float64}}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), (r["λ"], r["flip_fraction"]))
    end

    # Generate a color gradient from light blue to dark blue
    n_lines = length(grouped)
    blues = cgrad(:blues, n_lines, rev=true)  # from light to dark blue

    p = plot(xlabel="λ", ylabel="Fluctuant Fraction")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        sorted_data = sort(data, by=x->x[1])  # sort by λ
        lambdas = [x[1] for x in sorted_data]
        fluctuant_vals = [x[2] for x in sorted_data]
        plot!(p, lambdas, fluctuant_vals, label="N=$N, μ=$μ, k_avg=$k_avg", color=blues[i], linewidth=1.5, alpha=0.8)
        i += 1
    end
    plot!(p, legend=false)
    savefig(p, path)
end
#--Only for Structure Analysis------------------------------------------------------------------------------------------------------    
function plot_fluctuant_vs_trust_level_different_net(results, path)
    """Trust effect on system stability across structures"""
    
    # Group by network configuration
    grouped = Dict{Tuple{Int, Float64, Int}, Vector{Tuple{Float64, Float64}}}()
    for r in results
        k = (r["N"], r["μ"], r["k_avg"])
        push!(get!(grouped, k, []), (r["trust_level"], r["flip_fraction"]))
    end

    # Generate a color gradient from light red to dark red
    n_lines = length(grouped)
    reds = cgrad(:reds, n_lines, rev=true)  # from light to dark red

    p = plot(xlabel="Trust Level τ", ylabel="Fluctuant Fraction")
    i = 1
    for ((N, μ, k_avg), data) in grouped
        sorted_data = sort(data, by=x->x[1])  # sort by trust level
        trust_levels = [x[1] for x in sorted_data]
        fluctuant_vals = [x[2] for x in sorted_data]
        plot!(p, trust_levels, fluctuant_vals, label="N=$N, μ=$μ, k_avg=$k_avg", color=reds[i], linewidth=1.5, alpha=0.8)
        i += 1
    end
    plot!(p, legend=false)
    ylims!(p, (0, 1.05))
    savefig(p, path)
end
#-----------------------------------------------------------------------------------------------------------------------------------
end