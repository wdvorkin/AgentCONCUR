ns(l) = Int(n_s[l])
nr(l) = Int(n_s[l])
function remove_col_and_row(B,refbus)
    @assert size(B,1) == size(B,2)
    n = size(B,1)
    return B[1:n .!= refbus, 1:n .!= refbus]
end
function build_B̆(B̂inv,refbus)
    Nb = size(B̂inv,1)+1
    B̆ = zeros(Nb,Nb)
    for i in 1:Nb, j in 1:Nb
        if i < refbus && j < refbus
            B̆[i,j] = B̂inv[i,j]
        end
        if i > refbus && j > refbus
            B̆[i,j] = B̂inv[i-1,j-1]
        end
        if i > refbus && j < refbus
            B̆[i,j] = B̂inv[i-1,j]
        end
        if i < refbus && j > refbus
            B̆[i,j] = B̂inv[i,j-1]
        end
    end
    return B̆
end
function load_pwr_data()
    pwr_data = Dict()

    # load data sheets 
    node_data = CSV.read("NYISO_v2/Node.csv", DataFrame)
    gen_data =  CSV.read("NYISO_v2/Generator.csv", DataFrame)
    line_data = CSV.read("NYISO_v2/Line.csv", DataFrame)
    pwr_data[:node] = sort(Dict(node_data[:,:Node][i] => node_data[:,:name][i] for i in 1:11))

    pv_wd_generation = vcat(gen_data[gen_data[:,:gtype].=="Wind_LB",:gindex],
                            gen_data[gen_data[:,:gtype].=="Wind_OS",:gindex],
                            gen_data[gen_data[:,:gtype].=="PV",:gindex])
    conv_generation = setdiff(1:size(gen_data,1),pv_wd_generation)

    # system dims 
    D = size(node_data,1)
    N = size(node_data,1)
    G = size(conv_generation,1)
    E = size(line_data,1)
    pwr_data[:dims] = Dict(:D => D, :N => N, :G => G, :E => E)

    # demand data 
    d = node_data[!,:Pd]
    pwr_data[:load] = Dict(:d => d)

    # generation data 
    c = gen_data[conv_generation,14]
    C = diagm(gen_data[conv_generation,14]*0.00001)
    p̅ = gen_data[conv_generation,:Pmax]
    p̲ = gen_data[conv_generation,:Pmin]
    M_p = zeros(N,G)
    for i in 1:G
        M_p[gen_data[conv_generation,:location][i],i] = 1
    end
    pwr_data[:gen] = Dict(:c => c, :C => C, :p̅ => p̅, :p̲ => p̲, :M_p => M_p)

    # transmission data
    n_s = line_data[!,:from]
    n_r = line_data[!,:to]
    x = line_data[!,:x]
    f̅ = line_data[!,:s]
    ref = 1
    # Compute PTDF matrix
    B_line = zeros(E,N); B̃_bus = zeros(N,N); B = zeros(N,N)
    for n in 1:N
        for l in 1:E
            if n_s[l] == n
                B[n,n] += x[l]
                B_line[l,n] = x[l]
            end
            if n_r[l] == n
                B[n,n] += x[l]
                B_line[l,n] = -x[l]
            end
        end
    end
    for l in 1:E
        B[Int(n_s[l]),Int(n_r[l])] = - x[l]
        B[Int(n_r[l]),Int(n_s[l])] = - x[l]
    end
    B̃_bus = remove_col_and_row(B,ref)
    B̃_bus = inv(B̃_bus)
    B̃_bus = build_B̆(B̃_bus,ref)
    PTDF = B_line*B̃_bus
    pwr_data[:trans] = Dict(:n_s => n_s, :n_r => n_r, :x => x, :f̅ => f̅, :ref => ref, :PTDF => PTDF)

    return pwr_data
end
function data_split(dataset,set)
    N_train = set[:N_train]; 
    N_records = length(dataset);
    N_test = N_records - N_train;
    Random.seed!(256); train_keys = sample(1:N_records, N_train, replace = false); Random.seed!();
    test_keys = setdiff(1:N_records,train_keys);
    dataset_train = deepcopy(dataset);
    dataset_test = deepcopy(dataset);
    for i in test_keys
        delete!(dataset_train,i)
    end
    for i in train_keys
        delete!(dataset_test,i)
    end
    dataset_train = sort(Dict(keys(1:N_train) .=> values(dataset_train))); 
    dataset_test  = sort(Dict(keys(1:N_test)  .=> values(dataset_test)));
    return dataset_train, dataset_test
end
function performance_check(dataset,sol)
    N = length(dataset)
    net_dc_infeas = zeros(N)
    pwr_sy_infeas = zeros(N)
    latency_viol_status = zeros(N)
    latency_loss = zeros(N)
    opf_cost = zeros(N)
    for n in 1:N
        # compute the task shift
        task_shift = task_shift_from_reg_weights(sol,dataset[n])
        # check if feasible for a net of DCs 
        latency = latency_violation(net_dc_data,sol_lat_opt,task_shift[:φ])
        latency_loss[n] = latency[:latency_loss]
        if latency[:violation] <= -1
            latency_viol_status[n] = 1
        end
        # check if feasible for a power system 
        if latency_viol_status[n] != 1
            sol_opf = optimal_power_flow(pwr_data,dataset[n][:d],dataset[n][:r],task_shift)
        else
            sol_opf = optimal_power_flow(pwr_data,dataset[n][:d],dataset[n][:r],sol_lat_opt)
        end
        opf_cost[n] = sol_opf[:obj]

        # sol_opf[:status] == "INFEASIBLE" ? pwr_sy_infeas[n] = 1 : NaN
        # push!(load_shedding,sum(sol_opf[:l]))
    end
    return  Dict(:lat_viol => sum(latency_viol_status)/N*100, :lat_loss => mean(latency_loss), :pwr_sys_inf => sum(pwr_sy_infeas)/N*100, 
                 :opf_cost => mean(opf_cost))
end



function statistics(data)
    f = zeros(size(data[1][:f],1),length(data))
    d = zeros(length(data))
    π = zeros(size(data[1][:π],1),length(data))
    r = zeros(length(data))

    for i in 1:length(data)
        f[:,i] = abs.(data[i][:f])
        d[i] = sum(data[i][:d])
        π[:,i] = abs.(data[i][:π])
        r[i] = sum(data[i][:r])
    end
    f, d, π, r = vec(f), vec(d), vec(π), vec(r)
    @show round.(minimum(f),digits=2), round.(mean(f),digits=2), round.(maximum(f),digits=2)
    @show round.(minimum(d),digits=2), round.(mean(d),digits=2), round.(maximum(d),digits=2)
    @show round.(minimum(π),digits=2), round.(mean(π),digits=2), round.(maximum(π),digits=2)
    @show round.(minimum(r),digits=2), round.(mean(r),digits=2), round.(maximum(r),digits=2)
end





