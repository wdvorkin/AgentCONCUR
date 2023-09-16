# AgentCONCUR

Welcome to the online repository associated with the paper titled [*Agent Coordination via Contextual Regression (AgentCONCUR) for Data Center Flexibility*](https) submitted to the XXIII Power System Computational Conference. This repository is intended to provide supplementary materials, code, data, and resources to facilitate the understanding and replication of our results. If you use these materials, please be sure to cite the paper accordingly.

# Abstract 

A network of spatially distributed data centers can provide operational flexibility to power systems by shifting computing tasks among electrically remote locations. However, harnessing this flexibility in real-time through the standard optimization techniques is challenged by the need for sensitive operational datasets and substantial computational resources. To alleviate the data and computational requirements, this paper introduces a coordination mechanism based on contextual regression. This mechanism, abbreviated as AgentCONCUR, associates cost-optimal task shifts with public and trusted contextual data (e.g., real-time prices) and uses regression on this data as a coordination policy. Notably, regression-based coordination does not learn the optimal coordination actions from a labeled dataset. Instead, it exploits the optimization structure of the coordination problem to ensure feasible and cost-effective actions. A NYISO-based study reveals large coordination gains and the optimal features for the successful regression-based coordination.

![pic](https://github.com/wdvorkin/AgentCONCUR/assets/31773955/6417c329-40cd-4d6b-95a0-458a98986643)

# Citation

To be added 

# Usage

All models are implemented in Julia Language v1.8 using [JuMP](https://github.com/jump-dev/JuMP.jl), a modeling language for mathematical optimization, and the commercial optimization solvers [Mosek](https://github.com/MOSEK/Mosek.jl) and [Gurobi](https://github.com/jump-dev/Gurobi.jl), which need to be licensed (free for academic use). Also, make sure to activate the project environment using `Project.toml` and `Manifest.toml`.

The main execution file is `main.jl`, which computes the Power-NetDC coordination costs for varying NetDC penetration levels and the maximum allowable latency loss. To execute the file, clone this repository, open the terminal, navigate to the repository with `cd`, and type `julia main.jl -a 0.25`, which returns the results for a maximum allowable latency loss of 25%. To see all available options, use `julia main.jl --info`.








