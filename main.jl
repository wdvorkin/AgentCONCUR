using CSV, DataFrames, ArgParse
using LinearAlgebra, Statistics, StatsBase, Random, Distributions
using JuMP, Gurobi, Mosek, MosekTools
const GRB_ENV = Gurobi.Env()

# parse arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--alpha", "-a"
            help = "latency sub-optimality bound"
            arg_type = Float64
            default = 1.0
        "--net_dc_pen", "-p"
            help = "amount of net_dc electric demand to the peak load"
            arg_type = Float64
            default = 0.10 
        "--n_train", "-n"
            help = "number of training samples"
            arg_type = Int64
            default = 100
    end
    return parse_args(s)
end
args = parse_commandline()

# ### experiment settings 
set = Dict(:net_dc_pen => args["net_dc_pen"], 
:α => args["alpha"], 
:N_train => args["n_train"], 
:ε => 10,
:ϱ => 0.000001)

include("01_aux_fun.jl")
include("02_opt_fun.jl")
include("03_time_series_pre_processing.jl")
include("04_data_center_inputs.jl")

pwr_data = load_pwr_data()
load_comm, load_fore, renew_pwr = load_peak_data()
cost_summary = DataFrame(Share=Any[],Ideal=Any[],AgentConcure=Any[],BaselineReg=Any[],NoCoordination=Any[])

for q in 0.05:0.025:0.35
    set[:net_dc_pen] = q
    @info("current penetration level is $(set[:net_dc_pen])")
    ################################################################### 
    ###### load power system, peak, and datacenter network's parameters 
    ################################################################### 
    global net_dc_data = load_data_center_data(pwr_data,set)

    ###################################################################
    ###### optimal (baseline) latency solution ########################
    ###################################################################
    global sol_lat_opt = latency_opt(net_dc_data)

    ###################################################################
    ###### create a dataset for training and testing ##################
    ###################################################################
    N_records = length(load_comm)
    dataset = Dict()
    for n in 1:N_records
        record = Dict()
        record[:d] = load_comm[n]
        record[:r] = renew_pwr[n]
        sol_opf = optimal_power_flow(pwr_data,net_dc_data,load_comm[n],renew_pwr[n],sol_lat_opt)
        record[:π] = sol_opf[:π]
        record[:f] = sol_opf[:f]
        sol_optimal_load_shift = bilevel_optimization_of_load_shift(pwr_data,net_dc_data,load_comm[n],renew_pwr[n],sol_lat_opt)
        record[:φ_opt] =  sol_optimal_load_shift[:φ]
        dataset[n] = record
    end
    dataset = sort(dataset)
    dataset_train, dataset_test = data_split(dataset,set)

    ###################################################################
    ###### baseline regression training ###############################
    ###################################################################
    sol_reg_baselin = baseline_regression_training(dataset_train)

    ###################################################################
    ###### bi-level regression training ###############################
    ###################################################################
    sol_reg_bilevel = bilevel_regression_training(pwr_data,net_dc_data,dataset_train,sol_lat_opt)

    ###################################################################
    ###### cost summary ###############################################
    ###################################################################
    # IDEAL COORDINATION 
    cost_ideal_coordin = zeros(N_records-set[:N_train])
    for n in 1:(N_records-set[:N_train])
        sol_optimal_load_shift = bilevel_optimization_of_load_shift(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],sol_lat_opt)
        cost_ideal_coordin[n] = sol_optimal_load_shift[:obj]
    end

    # NO COORDINATION 
    cost_no_coordinati = zeros(N_records-set[:N_train])
    for n in 1:(N_records-set[:N_train])
        sol_opf = optimal_power_flow(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],sol_lat_opt)
        cost_no_coordinati[n] = sol_opf[:obj]
    end

    # BASELINE REG COORDINATION 
    cost_baseline_regr = zeros(N_records-set[:N_train])
    for n in 1:(N_records-set[:N_train])
        task_shift = task_shift_from_reg_weights(sol_reg_baselin,dataset_test[n])
        # OPF feasibility check
        opf_feas = opf_feas_check(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],task_shift)
        # latency feasibility check
        lat_feas = lat_feas_check(net_dc_data,sol_lat_opt,task_shift[:φ])
        # opf cost computation 
        if opf_feas[:status] == "OPTIMAL" && lat_feas[:status] == "OPTIMAL"
            sol_opf = optimal_power_flow(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],task_shift)
            cost_baseline_regr[n] = sol_opf[:obj]
        else
            cost_baseline_regr[n] = cost_no_coordinati[n] 
        end
    end

    # AgentCONCUR COORDINATION 
    cost_agent_concure = zeros(N_records-set[:N_train])
    for n in 1:(N_records-set[:N_train])
        task_shift = task_shift_from_reg_weights(sol_reg_bilevel,dataset_test[n])
        # OPF feasibility check
        opf_feas = opf_feas_check(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],task_shift)
        # latency feasibility check
        lat_feas = lat_feas_check(net_dc_data,sol_lat_opt,task_shift[:φ])
        # opf cost computation 
        if opf_feas[:status] == "OPTIMAL" && lat_feas[:status] == "OPTIMAL"
            sol_opf = optimal_power_flow(pwr_data,net_dc_data,dataset_test[n][:d],dataset_test[n][:r],task_shift)
            cost_agent_concure[n] = sol_opf[:obj]
        else
            cost_agent_concure[n] = cost_no_coordinati[n] 
        end
    end

    # mean costs
    cost_ideal_coordin_mean = mean(cost_ideal_coordin)
    cost_baseline_regr_mean = mean(cost_baseline_regr)
    cost_agent_concure_mean = mean(cost_agent_concure)
    cost_no_coordinati_mean = mean(cost_no_coordinati)

    push!(cost_summary,[q,cost_ideal_coordin_mean,cost_agent_concure_mean,cost_baseline_regr_mean,cost_no_coordinati_mean])
    display(cost_summary)
    CSV.write("cost_summary_alpha_$(set[:α]).csv", cost_summary)
end
