include("functions.jl")
using Plots

function plot_K_time_K_time_min_entropy_vsruido(noise_signal_range, d_values, filename)
    colors = [:magenta3, :green2, :blue, :brown3, :cyan, :black, :red]
    linestyle_minentropy = :dash
    
    first_plot = true  

    for (i, d) in enumerate(d_values)
        # Gráfica de "Nuestro método"
        K_time_vec = K_time.(Ref(4), Ref(d), noise_signal_range)
        label = "Nuestro método d = $d"  
        if first_plot
            plot(noise_signal_range, K_time_vec, label=label, color=colors[i], legend=true, ylim=(0, 2500),
                xminorgrid=true, yminorgrid=true, xminorticks=5, yminorticks=5
            )
            first_plot = false
        else
            plot!(noise_signal_range, K_time_vec, label=label, color=colors[i])
        end
        
        # Gráfica de "Min-entropía"
        K_time_min_entropy_vec = K_time_min_entropy.(Ref(d), noise_signal_range)
        label_min_entropy = "Min-entropía d = $d"
        plot!(noise_signal_range, K_time_min_entropy_vec, label=label_min_entropy, color=colors[i], linestyle=linestyle_minentropy)
    end
    
    xlabel!("Ruido/Señal")
    ylabel!("Tasa de clave (bits/s)")
    title!("")
    savefig(filename)
end

function plot_K_space_K_space_min_entropy(noise_signal_range, d_values, filename)
    colors = [:magenta3, :green2, :blue, :brown3, :blue, :cyan, :black, :red]
    linestyle_min_entropy = :dash
    
    first_plot = true  

    for (i, d) in enumerate(d_values)
        # Gráfica de "Nuestro método"
        K_space_vec = K_space.(Ref(4), Ref(d), noise_signal_range)
        label = "Nuestro método d = $d"  
        if first_plot
            plot(noise_signal_range, K_space_vec, label=label, color=colors[i], legend=true, ylim=(0, 2500),
                xminorgrid=true, yminorgrid=true, xminorticks=5, yminorticks=5
            )
            first_plot = false
        else
            plot!(noise_signal_range, K_space_vec, label=label, color=colors[i])
        end
        
        # Gráfica de "Min-entropía"
        K_space_min_entropy_vec = K_space_min_entropy.(Ref(d), noise_signal_range)
        label_min_entropy = "Min-entropía d = $d"
        plot!(noise_signal_range, K_space_min_entropy_vec, label=label_min_entropy, color=colors[i], linestyle=linestyle_min_entropy)
    end
    
    xlabel!("Ruido/Señal")
    ylabel!("Tasa de clave (bits/s)")
    title!("")
    savefig(filename)
end
