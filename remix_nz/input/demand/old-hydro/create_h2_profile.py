

import csv

def create_demand_profile(yearly_demands, names, list_of_percentage):
    # Create and write to the CSV file
    with open('H2_demand_profile_LUT.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Write the header
        header = ['Name', 'Regional Demand'] + [f'Hourly Demand {i+1}' for i in range(8760)]
        writer.writerow(header)
        
        # Iterate over each yearly demand and corresponding name
        for yearly_demand, name in zip(yearly_demands, names):
            # Calculate regional demands
            regional_demands = [yearly_demand * (percentage) for percentage in list_of_percentage]
            
            # Print the regional demands
            print(f"Regional Demands for {name}:", regional_demands)
            
            # Calculate hourly demands
            hourly_demands = [demand / 8760 for demand in regional_demands]
            
            # Print the hourly demands
            print(f"Hourly Demands for {name}:", hourly_demands)
            
            # Write the demand values for each region
            for regional_demand, hourly_demand in zip(regional_demands, hourly_demands):
                row = [name, regional_demand] + [hourly_demand] * 8760
                writer.writerow(row)




# demands from hadi
nodes_lst=["NIS","AKL","WTO","BOP", "HBY", "TRN", "CEN", "WEL", "NEL", "CAN", "OTG"]
shares_in_max_case=[0.0037, 0.2492, 0.0411, 0.0412, 0.0460, 0.2183, 0.0417, 0.1220, 0.0322,  0.1384, 0.0662]
shares_in_low_case=[0.0129, 0.0013, 0.0576, 0.0079, 0.0256, 0.7176, 0.0267, 0.0078, 0.0037, 0.0817, 0.0572]
list_of_percentage=[shares_in_max_case,shares_in_low_case]

#2050 from hadi
demand_kg = [860755, 243339] # max case  and low case 
demand_GWh= [demand * 33.33 for demand in demand_kg]
names_2050 =  ["Max case 2050", "Low case 2050"]

#from LUT
demand_lut = [42960, 88440, 105780]#H
names_lut =  ["Domestic 2030","Domestic 2050", "Exports 2050"]

#from leon / explained in his ppt
demand_pypsa = [654.557 , 1706.975, 20037.406] # 2030, 2050, 2050-cap(where there are limits in co2)
names_pypsa =  ["2030", "2050", "2050-cap"]
demand_pypsa_exports = [654.557, 20654.557, 40654.557, 200654.557] # just local demand, 20twh export, 40twh export, 200 twh export
names_pypsa_exports =  ["Domestic 2030", "20twh-export 2030", "40twh-export 2030", "200twh-export 2030"]


yearly_demands=[demand_GWh, demand_lut, demand_pypsa, demand_pypsa_exports] 
names=[names_2050, names_lut, names_pypsa, names_pypsa_exports] 
#create_demand_profile(yearly_demands[2], names[2], list_of_percentage[1])
create_demand_profile(yearly_demands[1], names[1], list_of_percentage[1])







