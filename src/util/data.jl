""" 
For data manipulations
"""

function mn_data(data, load_data, pv_data, inst_data; time_steps=8760)
    mn_model = PowerModels.replicate(data,time_steps,global_keys=Set{String}())
    power_base=data["baseMVA"]
    #mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in data["load"]
        mean_power = device["yearlyNetConsumption"]*365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        # print(device["source_id"])
        # print(device["source_id"] in names(load_data))
        if device["source_id"] in names(load_data)
            load_profile =  load_data[!,device["source_id"]]  #Pick the profile from csv
            load_profile = load_profile/1000 #scale from kWh to MW
            load_profile = load_profile/power_base #convert to per-uits
            pv_inc=inst_data[(inst_data[!,"Meter_name"].==device["source_id"]),:][!,"roof-top\norientation"]
            
            pv_prof = Matrix(pv_data[!,pv_inc])
            pv_gen=pv_prof.*inst_data[(inst_data[!,"Meter_name"].==device["source_id"]),:][!,"number of PV panels"]
            # print(pv_gen)
            pv_gen= pv_gen./(1e6*power_base)

            for step in 1:time_steps
                pd = load_profile[step]
                qd = pd/20
                #if length(device["phases"]) == 3   #Three phase connection
                    mn_model["nw"]["$(step)"]["load"][id_s]["pd"] = pd 
                    mn_model["nw"]["$(step)"]["load"][id_s]["qd"] = qd 
                
            end;


        
            #[mn_model["nw"]["$j"]["load"]["$i"]["pd"]=mn_model["nw"]["$j"]["load"]["$i"]["pd"]/12  for j=st_time[ind]:time_steps-8]
            [mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]=(mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]-pv_gen[j])/3  for j=1:time_steps]

        end;
    end
    return mn_model
end


""" 
For data prep for deterministic OPF on selected time_opf
"""

function mn_data_opf(data, load_data, pv_data, inst_data,time_opf)
    [data["PV"][d]["p_max"]=0 for d in keys(data["PV"])]
    mn_model = PowerModels.replicate(data,length(time_opf),global_keys=Set{String}())
    power_base=data["baseMVA"]
    
    #mn_model["data_model"] = MATHEMATICAL
    for (id_s,device) in data["load"]
        mean_power = device["yearlyNetConsumption"]*365/24 #Assuming yearlyNetConsumption contains total consumption of 365 days in kWh
        # print(device["source_id"])
        print(device["source_id"] in names(load_data))
        if device["source_id"] in names(load_data)
            load_profile =  load_data[!,device["source_id"]]  #Pick the profile from csv
            load_profile= coalesce.(load_profile, 0.0)
            load_profile = load_profile/1000 #scale from kWh to MW
            load_profile = load_profile/power_base #convert to per-uits
            pv_inc=inst_data[(inst_data[!,"Meter_name"].==device["source_id"]),:][!,"roof-top orientation"]
            
            pv_prof = Matrix(pv_data[!,pv_inc])
            pv_prof=coalesce.(pv_prof, 0.0)
            # pv=pv_prof.*inst_data[(inst_data[!,"Meter_name"].==device["source_id"]),:][!,"number of PV panels"]
            # print(pv_gen)
            pv_max=inst_data[(inst_data[!,"Meter_name"].==device["source_id"]),:][!,"number of PV panels"]
            pv_gen= pv_prof./(1e6*power_base)
            display(pv_max)
            
            for i in 1:length(time_opf)
                pd = load_profile[time_opf[i]]
                qd = pd/20
                #if length(device["phases"]) == 3   #Three phase connection
                    mn_model["nw"]["$(i)"]["load"][id_s]["pd"] = pd/3
                    mn_model["nw"]["$(i)"]["load"][id_s]["qd"] = qd/3
                    mn_model["nw"]["$(i)"]["PV"][id_s]["pd"] = pv_gen[time_opf[i]]/3*(1/0.22)
                    mn_model["nw"]["$(i)"]["PV"][id_s]["qd"] = 0
                    mn_model["nw"]["$(i)"]["PV"][id_s]["p_max"] = pv_max[1]*0.220
                    mn_model["nw"]["$(i)"]["PV"][id_s]["p_min"] = 0*0.220
                    mn_model["nw"]["$(i)"]["PV"][id_s]["p_no_max"] = pv_max[1]


                
            end;

        else
            for i in 1:length(time_opf)
                mn_model["nw"]["$(i)"]["PV"][id_s]["p_max"] = 0*0.220
                mn_model["nw"]["$(i)"]["PV"][id_s]["p_min"] = 0*0.220
                mn_model["nw"]["$(i)"]["PV"][id_s]["p_no_max"] = 0
                print("I am here")
            end
            # #[mn_model["nw"]["$j"]["load"]["$i"]["pd"]=mn_model["nw"]["$j"]["load"]["$i"]["pd"]/12  for j=st_time[ind]:time_steps-8]
            # [mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]=(mn_model["nw"]["$j"]["load"]["$id_s"]["pd"]-pv_gen[j])/3  for j=1:time_steps]

        end;
    end
    return mn_model
end
