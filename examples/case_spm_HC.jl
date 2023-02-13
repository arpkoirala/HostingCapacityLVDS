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
# This example is for Numerical Illustration C and D on Section IV of the paper      #
# In first Fig, OPF HC is compared with gPC-CC-OPF HC                          #
# Second figure is Comparision of 4 CC-OPF scenarios
• S1: The PV installations of consumers are in between 0 to 15 kWp. 
  The corresponding CC-OPF HC is 52.5 kWp.
• S2: The PV installations of consumers are equal, and in between 0 to 15 kWp.
  The corresponding CC-OPF HC of the test feeder is 44.8 kWp.
• S3: Consumer 1 and 2 in the test feeder can have PV installations
      between 0 and 10 kWp, consumer 3 has fixed PV installation of 5 kWp,
      and the remaining consumer can have PV installations between 0 and 15 kWp.
  The corresponding CC-OPF HC of the test feeder is 47.4 kWp.
• S4: All consumers have PV installations between 5 and 10 kWp.
     The corresponding CC-OPF HC is 47.4 kWp.
#Third fig is the tuning of moments based reformulation
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
ipopt_solver = Ipopt.Optimizer

# input
deg  = 2
aux  = true
red  = false

feeder ="Pola/1076069_1274129_mod_configuration.json" #feeder with 7 consumer for test
# data
file  = joinpath(BASE_DIR, "test/data/Spanish/")
load_file= "beta_lm_2016_8_6.csv"
pv_file = "beta_pm_2016_8_6.csv"

"""
# Code to plot the result of Fig. 5 without the scenario based method #
"""
data  = SPM.build_mathematical_model_single_phase(file, feeder,load_file, pv_file, t_s= 59)
result_hc= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red, stochastic=false)
s2 = Dict("output" => Dict("duals" => true))
result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
scatter([1,2,3,4,5,6,7],[result_hc_2["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="gPC-CC-OPF HC", figsize=(28,8))
scatter!([1,2,3,4,5,6,7],[result_hc["solution"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="OPF HC", figsize=(28,8))
plot!(xlabel="Device Id")
plot!(ylabel="PV size [kWp]")
plot!(title="Fig. 5: deterministic vs stochastic OPF")


"""
# Code to plot the result of Fig. 6 For different Scenario #
• S1: The PV installations of consumers are in between 0 to 15 kWp. The corresponding CC-OPF HC is 52.5 kWp.
• S2: The PV installations of consumers are equal, and in between 0 to 15 kWp. The corresponding CC-OPF HC of
      the test feeder is 44.8 kWp.
• S3: Consumer 1 and 2 in the test feeder can have PV installations between 0 and 10 kWp, consumer 3 has fixed
      PV installation of 5 kWp, and the remaining consumer can have PV installations between 0 and 15 kWp. The
      corresponding CC-OPF HC of the test feeder is 47.4 kWp.
• S4: All consumers have PV installations between 5 and 10 kWp. The corresponding CC-OPF HC is 47.4 kWp.
"""


result_hc_s1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
scatter([1,2,3,4,5,6,7],[result_hc_s1["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="S1", figsize=(28,8))

result_hc_s2= SPM.run_sopf_hc_equal_pv(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
scatter!([1,2,3,4,5,6,7],[result_hc_s2["objective"] for i=1:length(data["load"])],label="S2", figsize=(28,8))


[data["PV"]["$i"]["p_min"]=0 for i=1:2]    
[data["PV"]["$i"]["p_max"]=10 for i=1:2]
[data["PV"]["$i"]["p_min"]=5 for i=3:3]
[data["PV"]["$i"]["p_max"]=5 for i=3:3]

result_hc_s3= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
scatter!([1,2,3,4,5,6,7],[result_hc_s3["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="S3", figsize=(28,8))

[data["PV"]["$i"]["p_min"]=5 for i=1:7]    
[data["PV"]["$i"]["p_max"]=10 for i=1:7]
result_hc_s4= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
scatter!([1,2,3,4,5,6,7],[result_hc_s4["solution"]["nw"]["1"]["PV"]["$i"]["p_size"] for i=1:length(data["load"])],label="S4", figsize=(28,8))
plot!(xlabel="Device Id")
plot!(ylabel="PV size [kWp]")
plot!(title="Fig 6: Different Scenario")


"""
Fig. 8 Moment based reformulation

"""
data  = SPM.build_mathematical_model_single_phase(file, feeder,load_file, pv_file, t_s= 59)
result_hc_2= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)

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

    result_hc_1= SPM.run_sopf_hc(data, PM.IVRPowerModel, ipopt_solver, aux=aux, deg=deg, red=red; setting=s2)
end

scatter([0,1,2,3,4,5,6,7],[sum(sample(result_hc_2, "bus", i, "vs"; sample_size=100000).>(1.05^2))/100000 for i=1:length(data["bus"])], label="before tuning")
scatter!([0,1,2,3,4,5,6,7],[sum(sample(result_hc_1, "bus", i, "vs"; sample_size=100000).>1.05^2)/100000 for i=1:length(data["bus"])], label="after tuning")
plot!(xlabel="Bus Id")
plot!(ylabel="P(Ui>1.05)")
plot!(title="Fig 7: Probability of voltage in each node crossing 1.05 pu")
