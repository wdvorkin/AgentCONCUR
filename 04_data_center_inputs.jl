function load_data_center_data(pwr_data,set)
    N_dc = 5; M_dc = zeros(pwr_data[:dims][:N],N_dc);
    M_dc[1,1] = 1; M_dc[4,2] = 1; M_dc[3,3] = 1; M_dc[8,4] = 1; M_dc[10,5] = 1;
    n = N_dc; k = Int(n*(n-1)/2)
    A = zeros(n,k)
    A[1,1:4] .= 1
    for j in 1:4
        A[j+1,j] = -1
    end
    A[2,5:7] .= 1
    for j in 5:7
        A[j-2,j] = -1
    end
    A[3,8:9] .= 1
    for j in 8:9
        A[j-4,j] = -1
    end
    A[4,10] = 1
    for j in 10
        A[j-5,j] = -1
    end
    δ = round.(pwr_data[:load][:d] .* set[:net_dc_pen],digits=2)
    C = zeros(N_dc,pwr_data[:dims][:N])
    C[1,1] = 0;C[2,1] = 4;C[3,1] = 2;C[4,1] = 5;C[5,1] = 6;
    C[1,2] = 1;C[2,2] = 3;C[3,2] = 1;C[4,2] = 4;C[5,2] = 5;
    C[1,3] = 2;C[2,3] = 2;C[3,3] = 0;C[4,3] = 3;C[5,3] = 4;
    C[1,4] = 4;C[2,4] = 0;C[3,4] = 2;C[4,4] = 3;C[5,4] = 4;
    C[1,5] = 3;C[2,5] = 1;C[3,5] = 1;C[4,5] = 2;C[5,5] = 3;
    C[1,6] = 4;C[2,6] = 2;C[3,6] = 2;C[4,6] = 2;C[5,6] = 3;
    C[1,7] = 4;C[2,7] = 2;C[3,7] = 2;C[4,7] = 1;C[5,7] = 2;
    C[1,8] = 5;C[2,8] = 3;C[3,8] = 3;C[4,8] = 0;C[5,8] = 2;
    C[1,9] = 5;C[2,9] = 3;C[3,9] = 3;C[4,9] = 1;C[5,9] = 1;
    C[1,10] = 6;C[2,10] = 4;C[3,10] = 4;C[4,10] = 2;C[5,10] = 0;
    C[1,11] = 6;C[2,11] = 4;C[3,11] = 4;C[4,11] = 2;C[5,11] = 2;
    return dc_data =  Dict(:C => C, 
                    :δ => δ,
                    :A => A,
                    :M_dc => M_dc,
                    :N_dc => N_dc,
                    :pwr_bus => findall(x->x.==1,vec(sum(M_dc,dims=2))),
                    :k => k)
end