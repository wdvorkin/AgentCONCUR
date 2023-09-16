using CSV, DataFrames
using Plots
using LinearAlgebra

function load_peak_data()
    ### load commitment data 
    keys = readdir(pwd()*"/NYISO_v2/time_series/load_commitment")[2:end]
    load_commitment = Dict()
    for i in keys
        load_commitment[i] = readdir(pwd()*"/NYISO_v2/time_series/load_commitment/"*i)
    end

    n_keys = length(keys)
    n_file = [length(load_commitment[keys[i]]) for i in 1:n_keys]
    n_peak_hours = sum(n_file)

    peak_loads = []
    for i in 1:n_keys
        for j in 1:n_file[i]
            load = CSV.read(pwd()*"/NYISO_v2/time_series/load_commitment/"*keys[i]*"/"*load_commitment[keys[i]][j], DataFrame)
            time_stamp = unique(load[:,"Time Stamp"])[18]
            load = load[load[:,"Time Stamp"] .== time_stamp,:]
            load = Matrix(load[:,5:end])
            load = sum(load,dims=2)[:,1]'
            load = max.(0,load)
            push!(peak_loads,load)
        end
    end
    load_commitment = zeros(11,n_peak_hours)
    for i in 1:n_peak_hours
        load_commitment[:,i] = peak_loads[i]
    end
    # ["Capitl" "Centrl" "Dunwod" "Genese" "Hud Vl" "Longil" "Mhk Vl" "Millwd" "N.Y.C." "North" "West"]
    #     1         2        3        4         5       6        7        8       9        10      11
    load_commitment = load_commitment[vec([11 4 2 10 7 1 5 8 3 9 6]),:]

    ### load forecast data 
    keys = readdir(pwd()*"/NYISO_v2/time_series/load_forecast")[2:end]
    load_forecast = Dict()
    for i in keys
        load_forecast[i] = readdir(pwd()*"/NYISO_v2/time_series/load_forecast/"*i)
    end

    n_keys = length(keys)
    n_file = [length(load_forecast[keys[i]]) for i in 1:n_keys]
    n_peak_hours = sum(n_file)

    peak_loads = []
    for i in 1:n_keys
        for j in 1:n_file[i]
            load = CSV.read(pwd()*"/NYISO_v2/time_series/load_forecast/"*keys[i]*"/"*load_forecast[keys[i]][j], DataFrame)
            load = load[19,2:12]
            load = Vector(load)
            load = load'
            push!(peak_loads,load)
        end
    end
    load_forecast = zeros(11,n_peak_hours)
    for i in 1:n_peak_hours
        load_forecast[:,i] = peak_loads[i]
    end
    load_forecast = load_forecast[vec([11 4 2 10 7 1 5 8 3 9 6]),:]

    # plot(load_commitment[10,:],label="commitment")
    # plot!(load_forecast[10,:],label="forecast")

    ### renewables data
    pv_gen_2018 = CSV.read(pwd()*"/NYISO_v2/time_series/renewables/hourly_solar_gen_2018.csv", DataFrame)
    pv_gen_2019 = CSV.read(pwd()*"/NYISO_v2/time_series/renewables/hourly_solar_gen_2019.csv", DataFrame)
    pv_gen = vcat(pv_gen_2018, pv_gen_2019)
    pv_gen = vcat(pv_gen[pv_gen[:,:LOCAL_HOUR_END] .== 17,:],pv_gen[pv_gen[:,:LOCAL_HOUR_END] .== "17",:])
    pv_gen = pv_gen[:,:tot_solar_mwh]
    allocation = [15.8 15.8 15.8 15.8 15.8 0 0 0 0 0 20.9]./100
    pv_gen = Matrix((pv_gen .* allocation)')

    wp_gen_2018 = CSV.read(pwd()*"/NYISO_v2/time_series/renewables/hourly_wind_gen_2018.csv", DataFrame)
    wp_gen_2019 = CSV.read(pwd()*"/NYISO_v2/time_series/renewables/hourly_wind_gen_2019.csv", DataFrame)
    wp_gen = vcat(wp_gen_2018, wp_gen_2019)
    wp_gen = vcat(wp_gen[wp_gen[:,:local_hour_end] .== 17,:],wp_gen[wp_gen[:,:local_hour_end] .== "17",:])
    wp_gen = wp_gen[:,:tot_wind_mwh]
    allocation = [2.7 0 13.6 17.7 11.4 0 0 0 0 0 54.5]./100
    wp_gen = Matrix((wp_gen .* allocation)')

    renew_pwr = pv_gen .+ wp_gen

    renew_pwr_ = Dict()
    load_comm_ = Dict()
    load_fore_ = Dict()
    for i in 1:n_peak_hours
        renew_pwr_[i] = renew_pwr[:,i]
        load_comm_[i] = load_commitment[:,i]
        load_fore_[i] = load_forecast[:,i]
    end

    return sort(load_comm_), sort(load_fore_), sort(renew_pwr_)
end
