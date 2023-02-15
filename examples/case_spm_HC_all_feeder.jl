"""
################################################################################
#  A. Koirala, T. V. Aacker, M. U. Hashmi, R. D’hulst, and D. V. Hertem,       #
“Chance-constrained optimization based PV hosting capacity calculation using   #
general Polynomial Chaos,”                                                     #
Submitted to IEE Transactions of Power Systems.                                #
[Online]: https://tinyurl.com/yckbp58j                                         #
################################################################################
# Based on StochasticPowerModels.jl                                            #
# An extention package of PowerModels.jl for Stochastic (Optimal) Power Flow   #
# See http://github.com/timmyfaraday/StochasticPowerModels.jl                  #
################################################################################
# This example is for Case study on Section V of the paper      #
# The plot 1 shows the relationship between computation time and sources of uncertainties 
# in the feeder. Which is Fig. 11 in the paper above.
################################################################################

# variables
"""

using Pkg
Pkg.activate(".")
using JuMP
using Ipopt
using PowerModels
using StochasticPowerModels
using PowerModelsDistribution
using JSON
using DataFrames
using CSV
using Statistics
using Plots
# constants 
const PM = PowerModels
const SPM = StochasticPowerModels

# solvers
# ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
r = false
dir=BASE_DIR*"/test/data/Spanish/All_feeder"
all_feeder=CSV.read(dir*"/ts_all_feeder.csv",DataFrame,header=["conf", "ts","HC1","HC2"])
load_file= "beta_lm_2016_8_6.csv"
pv_file = "beta_pm_2016_8_6.csv"
all_feeder[!,"HC1"]=zeros(nrow(all_feeder))
all_feeder[!,"HC2"]=zeros(nrow(all_feeder))
all_feeder[!,"HC3"]=zeros(nrow(all_feeder))
hc1=[]
hc2=[]
hc3=[]
t_cc=[]
t_opf=[]
nodes=[]
consumers=[]
unc=[]
global i=0
for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    
    ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_cpu_time"=>500.0)
    s2 = Dict("output" => Dict("duals" => true))


    global i=i+1
    print("Feeder no: $i \n")
    #feeder="All_feeder/"*all_feeder[1,"conf"]
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder,load_file, pv_file, t_s= 59)
    push!(unc,length(data["sdata"]))
    push!(nodes,length(data["bus"]))
    push!(consumers,length(data["load"]))
    [data["bus"]["$i"]["vmin"]=0.8 for i=1:length(data["bus"])]
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2, stochastic=false)
   
    #start up
    [data["PV"]["$i"]["p_size_start"]=result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["PV"])]

    # ipopt_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>5, "max_iter"=>3000)
    result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2)
    e=1;
    m=1

    # if result_hc_2["termination_status"]== PM.LOCALLY_SOLVED
    # for i=1:length(data["bus"])
    #     if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
    #         l_old=data["bus"]["$i"]["λvmax"]
    #         m=sample(result_hc_2, "bus", i, "vs"; sample_size=100000)
    #         if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]^2
    #             data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.3
    #         elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]^2
    #                 data["bus"]["$i"]["λvmax"]=2*l_old-(data["bus"]["$i"]["vmax"]^2-mean(m))/std(m)+0.1
    #         elseif quantile(m,[0.95])[1]< 0.99*data["bus"]["$i"]["vmax"]^2
    #             data["bus"]["$i"]["λvmax"]= l_old-0.8
    #         elseif quantile(m,[0.95])[1]< 1*data["bus"]["$i"]["vmax"]^2
    #             data["bus"]["$i"]["λvmax"]= l_old-0.2
    #         end
    #     end
    # end
    if result_hc_2["termination_status"]== PM.LOCALLY_SOLVED
        for i=1:length(data["bus"])
            if -result_hc_2["solution"]["nw"]["1"]["bus"]["$i"]["dual_voltage_max"]>500
                l_old=data["bus"]["$i"]["λvmax"]
                m=sample(result_hc_2, "bus", i, "vs"; sample_size=100000)
                if quantile(m,[0.95])[1]>1.001*data["bus"]["$i"]["vmax"]^2
                    data["bus"]["$i"]["λvmax"]=l_old+0.4
                elseif quantile(m,[0.95])[1]>1.0001*data["bus"]["$i"]["vmax"]^2
                        data["bus"]["$i"]["λvmax"]=l_old+0.2
                elseif quantile(m,[0.95])[1]< 0.99*data["bus"]["$i"]["vmax"]^2
                    data["bus"]["$i"]["λvmax"]= l_old-0.6
                elseif quantile(m,[0.95])[1]< 1*data["bus"]["$i"]["vmax"]^2
                    data["bus"]["$i"]["λvmax"]= l_old-0.2
                end
            end
        end
    result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=r; setting=s2)

    
        push!(hc1,result_hc_2["objective"])
        push!(hc2,result_hc_1["objective"])
        push!(t_cc, result_hc_1["solve_time"])
    else
        push!(hc1,-1)
        push!(hc2,-1)
        push!(t_cc, result_hc_2["solve_time"])
    end
     
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"])
        push!(t_opf, result_hc["solve_time"])
    else
        push!(hc3,-1)
        push!(t_opf, result_hc["solve_time"])
    end
end

all_feeder[!,"HC0"]=hc1
all_feeder[!,"HC_CC"]=hc2
all_feeder[!,"HC_OPF"]=hc3
all_feeder[!,"t_opf"]=t_opf
all_feeder[!,"t_cc"]=t_cc
all_feeder[!,"consumers"]=consumers
all_feeder[!,"nodes"]= nodes
all_feeder[!, "unc"] = unc
CSV.write("PV_HC_feeders_with_start.csv",all_feeder)

scatter(all_feeder[all_feeder[!, :t_cc] .<500, :].unc, log10.(all_feeder[all_feeder[!, :t_cc] .<500, :].t_cc),label="gPC-CC-OPF", figsize=(28,8))
# scatter!([1,2,3,4,5,6,7],[result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="OPF HC", figsize=(28,8))
plot!(xlabel="Uncertainties [-]")
plot!(ylabel="log10 of computation time [sec]")
plot!(title="Fig. 11: Boxplot of log$_{10}$ of computational time for \text{gPC-\gls{cc}-\gls{opf}} in real LV feeders with respect to the number of uncertainties considered.")

"""
#deterministic

for b in eachrow(all_feeder)
    feeder="All_feeder/"*b.conf
    file  = joinpath(BASE_DIR, "test/data/Spanish/")
    data  = SPM.build_mathematical_model_single_phase(file, feeder, t_s= b.ts)

    s2 = Dict("output" => Dict("duals" => true))
    result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2, stochastic=false)
    e=1;
    if result_hc["termination_status"]== PM.LOCALLY_SOLVED
        push!(hc3,result_hc["objective"])
        #push!(hc2,result_hc_1["objective"])
    else
        push!(hc3,-1)
        #push!(hc2,-1)
    end
end

all_feeder[!,"HC3"]=hc3

"""