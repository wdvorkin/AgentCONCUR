function latency_opt(net_dc_data)
    (n,m) = size(net_dc_data[:C])
    model = Model(() -> Mosek.Optimizer())
    JuMP.set_silent(model)
    # model variables
    @variable(model, w[1:n,1:m] >= 0)
    @variable(model, ϑ[1:n] >= 0)
    # objective function
    @objective(model, Min, sum(net_dc_data[:C].*w) + 1/2*set[:ϱ]*sum(w.*w))
    # constraints
    @constraint(model, μ[j=1:m], sum(w[i,j] for i in 1:n) == net_dc_data[:δ][j])
    @constraint(model, ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    # solve model
    optimize!(model)
    # @info("Latency optimization terminates with status $(termination_status(model))")
    return Dict(:w => JuMP.value.(w), :ϑ => JuMP.value.(ϑ), :obj => JuMP.objective_value.(model), :lat => sum(net_dc_data[:C].*JuMP.value.(w)))
end

function load_shift_opt(dc_data,sol_lat,φ)
    (n,m) = size(dc_data[:C])
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    # model variables
    @variable(model, w[1:n,1:m])
    @variable(model, ϑ[1:n])
    @variable(model, t)
    # objective function
    @objective(model, Min, 1/2 * t)
    @constraint(model, [1/2;t;sum(dc_data[:C].*(w .- sol_lat[:w]))] in RotatedSecondOrderCone())
    # @constraint(model, [1/2;r;vec(w)] in RotatedSecondOrderCone())

    # constraints
    @constraint(model, μ[j=1:m], sum(w[i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    @constraint(model, κ, dc_data[:A] * φ .== ϑ .- sol_lat[:ϑ])
    @constraint(model, ω̲, w .>= 0)
    @constraint(model, γ, sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat[:w])) >= 0)
    # solve model
    optimize!(model)
    if "$(termination_status(model))" == "OPTIMAL"
        sol = Dict( :ϑ => JuMP.value.(ϑ), :w => JuMP.value.(w), :t => JuMP.value.(t),
                :obj => JuMP.objective_value.(model), :lat => sum(dc_data[:C].*JuMP.value.(w)),
                :status => "$(termination_status(model))")
        return sol
    else
        return Dict(:status => "$(termination_status(model))")
    end
end

function latency_violation(dc_data,sol_lat,φ)
    (n,m) = size(dc_data[:C])
    model = Model(() -> Mosek.Optimizer())
    JuMP.set_silent(model)
    # model variables
    @variable(model, w[1:n,1:m])
    @variable(model, ϑ[1:n])
    @variable(model, t)
    # objective function
    @objective(model, Min, 1/2 * t)
    @constraint(model, [1/2;t;sum(dc_data[:C].*(w .- sol_lat[:w]))] in RotatedSecondOrderCone())
    @constraint(model, μ[j=1:m], sum(w[i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    @constraint(model, κ, dc_data[:A] * φ .== ϑ .- sol_lat[:ϑ])
    @constraint(model, ω̲, w .>= 0)
    # @constraint(model, γ, sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat[:w])) >= 0)
    # solve model
    optimize!(model)

    viol = sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(JuMP.value.(w) .- sol_lat[:w]))

    return Dict(:violation => viol, :latency_loss => sum(dc_data[:C].*(JuMP.value.(w) .- sol_lat[:w])))
end

function optimal_power_flow(pwr_data,net_dc_data,load,renew_pwr,lat_opt_sol)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    @variable(model, p[1:pwr_data[:dims][:G]])
    @variable(model, l[1:pwr_data[:dims][:N]])

    @objective(model, Min, p'*pwr_data[:gen][:C]*p + pwr_data[:gen][:c]'*p + 1000*ones(pwr_data[:dims][:N])'*l)

    @constraint(model, λ_b,  ones(pwr_data[:dims][:N])'*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) >= 0)
    @constraint(model, λ_f̅,-pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, λ_f̲, pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, p .<= pwr_data[:gen][:p̅])
    @constraint(model, p .>= pwr_data[:gen][:p̲])
    @constraint(model, l .<= load .+ net_dc_data[:M_dc]*lat_opt_sol[:ϑ])
    @constraint(model, l .>= 0)

    optimize!(model)
    # @show "$(termination_status(model))"
    if "$(termination_status(model))" != "OPTIMAL"
        return Dict(:status => "$(termination_status(model))", :obj => 0)
    else
        return Dict(:status => "$(termination_status(model))", :p => JuMP.value.(p), :l => JuMP.value.(l), :obj => JuMP.objective_value(model), 
                :f =>  pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*JuMP.value.(p) .+ renew_pwr .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]),
                :π => JuMP.dual(λ_b)*ones(pwr_data[:dims][:N]) - pwr_data[:trans][:PTDF]' * (JuMP.dual.(λ_f̅) - JuMP.dual.(λ_f̲)))
    end
end



function opf_feas_check(pwr_data,net_dc_data,load,renew_pwr,lat_opt_sol)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    @variable(model, p[1:pwr_data[:dims][:G]])
    @variable(model, l[1:pwr_data[:dims][:N]] == 0)

    @objective(model, Min, sum(p) + sum(l))

    @constraint(model, λ_b,  ones(pwr_data[:dims][:N])'*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) >= 0)
    @constraint(model, λ_f̅,-pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, λ_f̲, pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*lat_opt_sol[:ϑ]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, p .<= pwr_data[:gen][:p̅])
    @constraint(model, p .>= pwr_data[:gen][:p̲])
    @constraint(model, l .<= load .+ net_dc_data[:M_dc]*lat_opt_sol[:ϑ])
    @constraint(model, l .>= 0)

    optimize!(model)
    return Dict(:status=> "$(termination_status(model))")
end

function lat_feas_check(dc_data,sol_lat,φ)
    (n,m) = size(dc_data[:C])
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    # model variables
    @variable(model, w[1:n,1:m])
    @variable(model, ϑ[1:n])

    # objective function
    @objective(model, Min, sum(w) + sum(ϑ))

    # constraints
    @constraint(model, μ[j=1:m], sum(w[i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    @constraint(model, κ, dc_data[:A] * φ .== ϑ .- sol_lat[:ϑ])
    @constraint(model, ω̲, w .>= 0)
    @constraint(model, γ, sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat[:w])) >= 0)
    # solve model
    optimize!(model)
    return Dict(:status => "$(termination_status(model))")
end

function bilevel_optimization_of_load_shift(pwr_data,dc_data,load,renew_pwr,sol_lat_opt)
    (n,m) = size(dc_data[:C])
    k = dc_data[:k]

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    @variable(model, φ[1:k])
    @variable(model, p[1:pwr_data[:dims][:G]])
    @variable(model, l[1:pwr_data[:dims][:N]])
    @variable(model, w[1:n,1:m]>=0)
    @variable(model, ϑ[1:n])
    @variable(model, μ[1:m])
    @variable(model, ν[1:n])
    @variable(model, κ[1:n])
    @variable(model, ω̲[1:n,1:m]>=0)
    @variable(model, γ>=0)

    # optimal power flow cost
    @objective(model, Min, p'*pwr_data[:gen][:C]*p + pwr_data[:gen][:c]'*p + 1000*ones(pwr_data[:dims][:N])'*l)

    # optimal power flow constraints
    @constraint(model, ones(pwr_data[:dims][:N])'*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) == 0)
    @constraint(model,     -pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) .>= -pwr_data[:trans][:f̅])
    @constraint(model,      pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) .>= -pwr_data[:trans][:f̅])
    @constraint(model, p .<= pwr_data[:gen][:p̅])
    @constraint(model, p .>= pwr_data[:gen][:p̲])
    @constraint(model, l .<= load .+ net_dc_data[:M_dc]*ϑ)
    @constraint(model, l .>= 0)

    # net_dc primal feasibility
    @constraint(model, con_dc_1[j=1:m], sum(w[i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, con_dc_2[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    @constraint(model, dc_data[:A] * φ .== ϑ .- sol_lat_opt[:ϑ])
    @constraint(model, sum(dc_data[:C].*sol_lat_opt[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat_opt[:w])) >= 0)

    # stat conditions
    for i in 1:n, j in 1:m
        @constraint(model, dc_data[:C][i,j]*sum(dc_data[:C].*(w .- sol_lat_opt[:w])) .- μ[j] .- ν[i] .- ω̲[i,j] .-γ*dc_data[:C][i,j] .== 0)
    end
    @constraint(model, ν .+ κ .== 0)

    # complementarity slackness
    for i in 1:n, j in 1:m
        @constraint(model, [w[i,j];ω̲[i,j]] in SOS1())
    end
    @constraint(model, [γ; sum(dc_data[:C].*sol_lat_opt[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat_opt[:w]))] in SOS1())

    # solve model
    optimize!(model)
    if "$(termination_status(model))" != "OPTIMAL"
        @warn("Bilevel load shift optimization termination status is not optimal")
    end

    return Dict(:φ => JuMP.value.(φ), :ϑ => JuMP.value.(ϑ), :obj => JuMP.objective_value(model), :status => termination_status(model))
end

function baseline_regression_training(dataset)
    n_f = size(dataset[1][:f],1)
    n_π = size(dataset[1][:π],1)
    n_d = size(dataset[1][:d],1)
    n_r = size(dataset[1][:r],1)


    (n,m) = size(net_dc_data[:C])
    k = Int(n*(n-1)/2)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)

    @variable(model, β[1:(k+k*n_f+k*n_π+k*n_d+k*n_r)])
    @variable(model, v[1:(k+k*n_f+k*n_π+k*n_d+k*n_r)])
    @variable(model, β_0[1:k])
    @variable(model, β_f[1:k,1:n_f])
    @variable(model, β_π[1:k,1:n_π])
    @variable(model, β_d[1:k,1:n_d])
    @variable(model, β_r[1:k,1:n_r])
    # regularization variables 
    @variable(model, t[1:length(dataset)])
    @objective(model, Min, sum(t))
    @constraint(model, con[n=1:length(dataset)], [t[n];β_0 .+ β_f*dataset[n][:f] .+ β_π*dataset[n][:π] .+ β_d*dataset[n][:d] .+ β_r*dataset[n][:r] .- dataset[n][:φ_opt]] in SecondOrderCone())

    # l1-regularization
    @constraint(model, β .== vcat(vec(β_0),vec(β_f),vec(β_π),vec(β_d),vec(β_r)))
    @constraint(model, v .>= -β)
    @constraint(model, v .>=  β)
    @constraint(model, sum(v) <= set[:ε])

    optimize!(model)
    if "$(termination_status(model))" != "OPTIMAL"
        @warn("Baseline regression training termination status is not optimal")
    end


    return Dict(:β => JuMP.value.(β), :β_0 => JuMP.value.(β_0), :β_f => JuMP.value.(β_f), :β_π => JuMP.value.(β_π), :β_d => JuMP.value.(β_d), :β_r => JuMP.value.(β_r))
end

function task_shift_from_reg_weights(sol,data)
    φ̂ = sol[:β_0] .+ sol[:β_f]*data[:f] .+ sol[:β_π]*data[:π] .+ sol[:β_d]*data[:d] .+ sol[:β_r]*data[:r]
    return Dict(:ϑ => net_dc_data[:A]*φ̂ .+ sol_lat_opt[:ϑ], :φ => φ̂)
end




function load_shift_opt_relaxed(dc_data,sol_lat,φ)
    (n,m) = size(dc_data[:C])
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    JuMP.set_silent(model)
    # model variables
    @variable(model, w[1:n,1:m])
    @variable(model, ϑ[1:n])
    @variable(model, t)
    # objective function
    @objective(model, Min, 1/2 * t)
    @constraint(model, [1/2;t;sum(dc_data[:C].*(w .- sol_lat[:w]))] in RotatedSecondOrderCone())
    # @constraint(model, [1/2;r;vec(w)] in RotatedSecondOrderCone())

    # constraints
    @constraint(model, μ[j=1:m], sum(w[i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
    @constraint(model, κ, dc_data[:A] * φ .== ϑ .- sol_lat[:ϑ])
    @constraint(model, ω̲, w .>= 0)
    # @constraint(model, γ, sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(w .- sol_lat[:w])) >= 0)
    # solve model
    optimize!(model)

    if "$(termination_status(model))" == "OPTIMAL"
        lat_viol_flag = 0
        if sum(dc_data[:C].*sol_lat[:w])*set[:α]-sum(dc_data[:C].*(JuMP.value.(w) .- sol_lat[:w])) >= 0.001
            lat_viol_flag = 1
        end
        sol = Dict(:latency => lat_viol_flag)
    else
        sol = Dict(:latency => 1)
    end
    return sol
end



function bilevel_regression_training(pwr_data,dc_data,dataset,sol_lat_opt)
    N = length(dataset)
    (n,m) = size(dc_data[:C])
    k = dc_data[:k]
    n_f = size(dataset[1][:f],1)
    n_π = size(dataset[1][:π],1)
    n_d = size(dataset[1][:d],1)
    n_r = size(dataset[1][:r],1)

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    # JuMP.set_silent(model)
    @variable(model, φ[1:N,1:k])
    @variable(model, p[1:N,1:pwr_data[:dims][:G]])
    @variable(model, l[1:N,1:pwr_data[:dims][:N]])
    @variable(model, w[1:N,1:n,1:m]>=0)
    @variable(model, ϑ[1:N,1:n])
    @variable(model, μ[1:N,1:m])
    @variable(model, ν[1:N,1:n])
    @variable(model, κ[1:N,1:n])
    @variable(model, ω̲[1:N,1:n,1:m]>=0)
    @variable(model, γ[1:N]>=0)

    @variable(model, β[1:(k+k*n_f+k*n_π+k*n_d+k*n_r)])
    @variable(model, v[1:(k+k*n_f+k*n_π+k*n_d+k*n_r)])
    @variable(model, β_0[1:k])
    @variable(model, β_f[1:k,1:n_f])
    @variable(model, β_π[1:k,1:n_π])
    @variable(model, β_d[1:k,1:n_d])
    @variable(model, β_r[1:k,1:n_r])

    # optimal power flow cost
    @objective(model, Min, 1/N*sum(p[s,:]'*pwr_data[:gen][:C]*p[s,:] + pwr_data[:gen][:c]'*p[s,:] + 1000*ones(pwr_data[:dims][:N])'*l[s,:] for s in 1:N))

    # regression constraints 
    @constraint(model, con_reg[s=1:N], φ[s,:] .== β_0 .+ β_f*dataset[s][:f] .+ β_d*dataset[s][:d] .+ β_π*dataset[s][:π] .+ β_r*dataset[s][:r])
    @constraint(model, β .== vcat(vec(β_0),vec(β_f),vec(β_π),vec(β_d),vec(β_r)))
    @constraint(model, v .>= -β)
    @constraint(model, v .>=  β)
    @constraint(model, sum(v) <= set[:ε])

    # optimal power flow constraints
    @constraint(model, con_1[s=1:N], ones(pwr_data[:dims][:N])'*  (pwr_data[:gen][:M_p]*p[s,:] .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) == 0)
    @constraint(model, con_2[s=1:N],     -pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p[s,:] .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, con_3[s=1:N],      pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p[s,:] .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) .>= -pwr_data[:trans][:f̅])
    @constraint(model, con_4[s=1:N], p[s,:] .<= pwr_data[:gen][:p̅])
    @constraint(model, con_5[s=1:N], p[s,:] .>= pwr_data[:gen][:p̲])
    @constraint(model, con_6[s=1:N], l[s,:] .<= dataset[s][:d] .+ net_dc_data[:M_dc]*ϑ[s,:])
    @constraint(model, con_7[s=1:N], l[s,:] .>= 0)

    # net_dc primal feasibility
    @constraint(model, con_dc_1[s=1:N,j=1:m], sum(w[s,i,j] for i in 1:n) == dc_data[:δ][j])
    @constraint(model, con_dc_2[s=1:N,i=1:n], sum(w[s,i,j] for j in 1:m) == ϑ[s,i])
    @constraint(model, con_dc_3[s=1:N], dc_data[:A] * φ[s,:] .== ϑ[s,:] .- sol_lat_opt[:ϑ])
    @constraint(model, con_dc_4[s=1:N], sum(dc_data[:C].*sol_lat_opt[:w])*set[:α]-sum(dc_data[:C].*(w[s,:,:] .- sol_lat_opt[:w])) >= 0)

    # stat conditions
    for s in 1:N, i in 1:n, j in 1:m
        @constraint(model, dc_data[:C][i,j]*sum(dc_data[:C].*(w[s,:,:] .- sol_lat_opt[:w])) .- μ[s,j] .- ν[s,i] .- ω̲[s,i,j] .-γ[s]*dc_data[:C][i,j] .== 0)
    end
    @constraint(model, ν .+ κ .== 0)

    # complementarity slackness
    for s in 1:N, i in 1:n, j in 1:m
        @constraint(model, [w[s,i,j];ω̲[s,i,j]] in SOS1())
    end
    @constraint(model, con_[s=1:N], [γ[s]; sum(dc_data[:C].*sol_lat_opt[:w])*set[:α]-sum(dc_data[:C].*(w[s,:,:] .- sol_lat_opt[:w]))] in SOS1())

    # solve model
    optimize!(model)
    if "$(termination_status(model))" != "OPTIMAL"
        @warn("Bilevel load shift optimization termination status is not optimal")
    end

    return Dict(:φ => JuMP.value.(φ), :ϑ => JuMP.value.(ϑ), :obj => JuMP.objective_value(model), :status => termination_status(model),
                :β => JuMP.value.(β), :β_0 => JuMP.value.(β_0), :β_f => JuMP.value.(β_f), :β_π => JuMP.value.(β_π), :β_d => JuMP.value.(β_d), :β_r => JuMP.value.(β_r),
                :CPU_time => JuMP.solve_time(model))
end


# function complementarity_check(task_shift,net_dc_data)
#     return sum(abs.(task_shift[:w].*task_shift[:ω̲])) + abs.(task_shift[:γ]*(sum(net_dc_data[:C].*lat_opt_sol[:w])*set[:α]-sum(net_dc_data[:C].*(task_shift[:w] .- lat_opt_sol[:w]))))
# end

# function bilevel_task_shift_optimization(pwr_data,load,renew_pwr,lat_opt_sol)
#     (n,m) = size(net_dc_data[:C])
#     model = Model(() -> Gurobi.Optimizer(GRB_ENV))
#     JuMP.set_silent(model)
#     # coordination variables 
#     @variable(model, φ[1:Int(n*(n-1)/2)])
#     # OPF variables
#     @variable(model, p[1:pwr_data[:dims][:G]])
#     @variable(model, l[1:pwr_data[:dims][:N]])
#     # NetDC primal variables 
#     @variable(model, w[1:n,1:m])
#     @variable(model, ϑ[1:n])
#     # NetDC dual variables 
#     @variable(model, μ[1:m])
#     @variable(model, ν[1:n])
#     @variable(model, κ[1:n])
#     @variable(model, ω̲[1:n,1:m])
#     @variable(model, γ)
#     # # regularization variables 
#     # @variable(model, t[1:Int(n*(n-1)/2)])

#     @objective(model, Min,  p'*pwr_data[:gen][:C]*p + pwr_data[:gen][:c]'*p 
#                                 + 1000*ones(pwr_data[:dims][:N])'*l)
#     # @constraint(model, φ .== 0)
#     # OPF constraints 
#     @constraint(model, λ_b,  ones(pwr_data[:dims][:N])'*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) >= 0)
#     @constraint(model, λ_f̅,    -pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) .>= -pwr_data[:trans][:f̅])
#     @constraint(model, λ_f̲,     pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p .+ renew_pwr .+ l .- load .- net_dc_data[:M_dc]*ϑ) .>= -pwr_data[:trans][:f̅])
#     @constraint(model, p .<= pwr_data[:gen][:p̅])
#     @constraint(model, p .>= pwr_data[:gen][:p̲])
#     @constraint(model, l .<= load)
#     @constraint(model, l .>= 0)

#     # NetDC constraints 
#     # primal feasibility 
#     @constraint(model, con_μ[j=1:m], sum(w[i,j] for i in 1:n) == net_dc_data[:δ][j])
#     @constraint(model, con_ν[i=1:n], sum(w[i,j] for j in 1:m) == ϑ[i])
#     @constraint(model, con_κ, net_dc_data[:A] * φ .== ϑ .- lat_opt_sol[:ϑ])
#     @constraint(model, con_ω̲, w .>= 0)
#     @constraint(model, con_γ, sum(net_dc_data[:C].*(w .- lat_opt_sol[:w])) <= sum(net_dc_data[:C].*lat_opt_sol[:w])*set[:α])
#     # dual feasibility
#     @constraint(model, ω̲ .>= 0)
#     @constraint(model, γ  >= 0)
#     # stationarity conditions
#     @constraint(model, ν .+ κ .== 0)
#     @constraint(model, ∂w[i=1:n,j=1:m], net_dc_data[:C][i,j] * (sum(net_dc_data[:C].*(w .- lat_opt_sol[:w]))) - μ[j] - ν[i] - ω̲[i,j] - γ * net_dc_data[:C][i,j]  .== 0)
#     # complementarity slackness
#     for i in 1:n, j in 1:m 
#         @constraint(model, [w[i,j];ω̲[i,j]] in SOS1())
#     end
#     @constraint(model, [γ;sum(net_dc_data[:C].*lat_opt_sol[:w])*set[:α]-sum(net_dc_data[:C].*(w .- lat_opt_sol[:w]))] in SOS1())
#     optimize!(model)
#     @show termination_status(model)

#     return Dict(
#         :p => JuMP.value.(p), 
#         :l => JuMP.value.(l),
#         :w => JuMP.value.(w), 
#         :ϑ => JuMP.value.(ϑ), 
#         :μ => JuMP.value.(μ), 
#         :ν => JuMP.value.(ν), 
#         :κ => JuMP.value.(κ), 
#         :ω̲ => JuMP.value.(ω̲), 
#         :γ => JuMP.value.(γ), 
#         :φ => JuMP.value.(φ), 
#         :obj => JuMP.objective_value(model))
# end

# function compute_load_shift(sol,data)
#     return sol[:β_0] .+ sol[:β_f] * data[:f] .+ sol[:β_π] * data[:π] .+ sol[:β_d] * data[:d] .+ sol[:β_r] * data[:r]
# end

# function base_line_regression(dataset)
#     n_f = size(dataset[1][:f],1)
#     n_π = size(dataset[1][:π],1)
#     n_d = size(dataset[1][:d],1)
#     n_r = size(dataset[1][:r],1)
#     k = length(dataset)
#     (n,m) = size(net_dc_data[:C])
#     model = Model(() -> Gurobi.Optimizer(GRB_ENV))
#     JuMP.set_silent(model)
#     o = Int(n*(n-1)/2)
#     @variable(model, β[1:(o+o*n_f+o*n_π+o*n_d+o*n_r)])
#     @variable(model, v[1:(o+o*n_f+o*n_π+o*n_d+o*n_r)])
#     @variable(model, β_0[1:o])
#     @variable(model, β_f[1:o,1:n_f])
#     @variable(model, β_π[1:o,1:n_π])
#     @variable(model, β_d[1:o,1:n_d])
#     @variable(model, β_r[1:o,1:n_r])
#     # regularization variables 
#     @variable(model, t[1:k])
#     @objective(model, Min, sum(t))
#     @constraint(model, con[n=1:k], [t[n];β_0 .+ β_f*dataset[n][:f] .+ β_π*dataset[n][:π] .+ β_d*dataset[n][:d] .+ β_r*dataset[n][:r] .- dataset[n][:φ]] in SecondOrderCone())

#     # l1-regularization
#     @constraint(model, β .== vcat(vec(β_0),vec(β_f),vec(β_π),vec(β_d),vec(β_r)))
#     @constraint(model, v .>= -β)
#     @constraint(model, v .>=  β)
#     @constraint(model, sum(v) <= set[:ε])

#     optimize!(model)

#     # @show termination_status(model)
#     return Dict(:obj => JuMP.objective_value(model), :β_0 => JuMP.value.(β_0), :β_f => JuMP.value.(β_f), :β_π => JuMP.value.(β_π), :β_d => JuMP.value.(β_d), :β_r => JuMP.value.(β_r), 
#                         :β => vcat(vec(JuMP.value.(β_0)),vec(JuMP.value.(β_f)),vec(JuMP.value.(β_π)),vec(JuMP.value.(β_d)),vec(JuMP.value.(β_r))))
# end

# function bilevel_regression_training(pwr_data,dataset,lat_opt_sol)
#     n_f = size(dataset[1][:f],1)
#     n_π = size(dataset[1][:π],1)
#     n_d = size(dataset[1][:d],1)
#     n_r = size(dataset[1][:r],1)
#     k = length(dataset)
#     (n,m) = size(net_dc_data[:C])
#     model = Model(() -> Gurobi.Optimizer(GRB_ENV))
#     JuMP.set_silent(model)

#     # coordination variables 
#     @variable(model, φ[1:k,1:Int(n*(n-1)/2)])
#     # OPF variables
#     @variable(model, p[1:k,1:pwr_data[:dims][:G]])
#     @variable(model, l[1:k,1:pwr_data[:dims][:N]])
#     # NetDC primal variables 
#     @variable(model, w[1:k,1:n,1:m])
#     @variable(model, ϑ[1:k,1:n])
#     # NetDC dual variables 
#     @variable(model, μ[1:k,1:m])
#     @variable(model, ν[1:k,1:n])
#     @variable(model, κ[1:k,1:n])
#     @variable(model, ω̲[1:k,1:n,1:m])
#     @variable(model, γ[1:k])
#     # regression variables 
#     o = Int(n*(n-1)/2)
#     @variable(model, β[1:(o+o*n_f+o*n_π+o*n_d+o*n_r)])
#     @variable(model, t[1:(o+o*n_f+o*n_π+o*n_d+o*n_r)])
#     @variable(model, β_0[1:o])
#     @variable(model, β_f[1:o,1:n_f])
#     @variable(model, β_π[1:o,1:n_π])
#     @variable(model, β_d[1:o,1:n_d])
#     @variable(model, β_r[1:o,1:n_r])

#     @objective(model, Min,  1/k*sum(p[s,:]'*pwr_data[:gen][:C]*p[s,:] + pwr_data[:gen][:c]'*p[s,:] 
#                                 + 1000*ones(pwr_data[:dims][:N])'*l[s,:] for s in 1:k))

#     # regression definition
#     @constraint(model, con_reg_1[s=1:k], φ[s,:] .== β_0 .+ β_f*dataset[s][:f] .+ β_π*dataset[s][:π] .+ β_d*dataset[s][:d] .+ β_r*dataset[s][:r])
#     @constraint(model, β .== vcat(vec(β_0),vec(β_f),vec(β_π),vec(β_d),vec(β_r)))

#     # @constraint(model,φ[:,5:end] .== 0)

#     # l1-regularization for feature selection 
#     @constraint(model, t .>= -β)
#     @constraint(model, t .>=  β)
#     @constraint(model, sum(t) <= set[:ε])

#     # OPF constraints 
#     @constraint(model, con_1[s=1:k],  ones(pwr_data[:dims][:N])'*(pwr_data[:gen][:M_p]*p[s,:]  .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) >= 0)
#     @constraint(model, con_2[s=1:k],    -pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p[s,:]  .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) .>= -pwr_data[:trans][:f̅])
#     @constraint(model, con_3[s=1:k],     pwr_data[:trans][:PTDF]*(pwr_data[:gen][:M_p]*p[s,:]  .+ dataset[s][:r] .+ l[s,:] .- dataset[s][:d] .- net_dc_data[:M_dc]*ϑ[s,:]) .>= -pwr_data[:trans][:f̅])
#     @constraint(model, con_4[s=1:k], p[s,:]  .<= pwr_data[:gen][:p̅])
#     @constraint(model, con_5[s=1:k], p[s,:]  .>= pwr_data[:gen][:p̲])
#     @constraint(model, con_6[s=1:k], l[s,:]  .<= dataset[s][:d])
#     @constraint(model, con_7[s=1:k], l[s,:]  .>= 0)

#     # NetDC constraints
#     # primal feasibility
#     @constraint(model, con_μ[s=1:k,j=1:m], sum(w[s,i,j] for i in 1:n) == net_dc_data[:δ][j])
#     @constraint(model, con_ν[s=1:k,i=1:n], sum(w[s,i,j] for j in 1:m) == ϑ[s,i])
#     @constraint(model, con_κ[s=1:k], net_dc_data[:A] * φ[s,:] .== ϑ[s,:] .- lat_opt_sol[:ϑ])
#     @constraint(model, con_ω̲[s=1:k], w[s,:,:] .>= 0)
#     @constraint(model, con_γ[s=1:k], sum(net_dc_data[:C].*(w[s,:,:] .- lat_opt_sol[:w])) <= sum(net_dc_data[:C].*lat_opt_sol[:w])*set[:α])
#     # dual feasibility
#     @constraint(model, con_df_1[s=1:k], ω̲[s,:,:] .>= 0)
#     @constraint(model, con_df_2[s=1:k], γ[s]  >= 0)
#     # stationarity conditions
#     @constraint(model, con_sc_1[s=1:k], ν[s,:] .+ κ[s,:] .== 0)
#     @constraint(model, con_sc_2[s=1:k,i=1:n,j=1:m], net_dc_data[:C][i,j] * (sum(net_dc_data[:C].*(w[s,:,:] .- lat_opt_sol[:w]))) - μ[s,j] - ν[s,i] - ω̲[s,i,j] - γ[s] * net_dc_data[:C][i,j]  .== 0)
#     # complementarity slackness
#     @constraint(model, con_cs_1[s=1:k,i=1:n,j=1:m], [w[s,i,j];ω̲[s,i,j]] in SOS1())
#     @constraint(model, con_cs_2[s=1:k], [γ[s];sum(net_dc_data[:C].*lat_opt_sol[:w])*set[:α]-sum(net_dc_data[:C].*(w[s,:,:] .- lat_opt_sol[:w]))] in SOS1())
#     optimize!(model)
#     @show termination_status(model)

#     return Dict(
#         :cost_per_scenario => [JuMP.value.(p[s,:])'*pwr_data[:gen][:C]*JuMP.value.(p[s,:]) + pwr_data[:gen][:c]'*JuMP.value.(p[s,:]) + 1000*ones(pwr_data[:dims][:N])'*JuMP.value.(l[s,:]) for s in 1:k],
#         :p => JuMP.value.(p), 
#         :l => JuMP.value.(l),
#         :w => JuMP.value.(w), 
#         :ϑ => JuMP.value.(ϑ), 
#         :μ => JuMP.value.(μ), 
#         :ν => JuMP.value.(ν), 
#         :κ => JuMP.value.(κ), 
#         :ω̲ => JuMP.value.(ω̲), 
#         :γ => JuMP.value.(γ), 
#         :φ => JuMP.value.(φ), 
#         :obj => JuMP.objective_value(model),
#         :β => JuMP.value.(β), 
#         :β_0 => JuMP.value.(β_0), 
#         :β_f => JuMP.value.(β_f), 
#         :β_π => JuMP.value.(β_π), 
#         :β_d => JuMP.value.(β_d), 
#         :β_r => JuMP.value.(β_r),
#         :t => JuMP.value.(t))
# end

# compute_new_net_dc_demand(φ) = Dict(:ϑ => net_dc_data[:A]*φ + lat_opt_sol[:ϑ])



