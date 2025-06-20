# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 08:24:18 2024
Modified from Oct-Dec 2024

@author: npw7gn, phw6bv
"""
import scipy
import numpy as np
import pandas as pd
import NetfluxODE_params as params
import NetfluxODE as file

[speciesNames, tau, ymax, y0, w, n, EC50] = params.loadParams()

#%%
tspan = [0, 40] 

sol_dict = {}
t_list = []
y_list = []

sol = scipy.integrate.solve_ivp(file.ODEfunc,tspan,y0,args=(tau,ymax,w,n,EC50,))
sol_dict[1] = sol
t_list.append(sol.t)
y_list.append(sol.y)

t = np.hstack(t_list)
y = np.hstack(y_list)

y_plot = y.T

ss = y_plot[-1,:] # steady state
#print(ss)

#%%
#create a for loop to analyze the change to the phenotypes and look for the highest change
#print(speciesNames[36])
ymax_knockdown = ymax.copy()
ymax_knockdown[19] = 0 #knowck down of PI3k, whihc is along the path ffor PI3KCA since all the drugs point to it


tspan = [0, 40] 

sol2_dict = {}
t2_list = []
y2_list = []

sol2 = scipy.integrate.solve_ivp(file.ODEfunc,tspan,y0,args=(tau,ymax_knockdown,w,n,EC50,))
sol2_dict[1] = sol2
t2_list.append(sol2.t)
y2_list.append(sol2.y)

t2 = np.hstack(t2_list)
y2 = np.hstack(y2_list)

y2_plot = y2.T

ss_pertubation = y2_plot[-1,:] # steady state after pertubation

#print(ss_pertubation)

change_in_activity = ss_pertubation - ss
# print(change_in_activity)

# %%
# Identify the most affected biomolecules by magnitude of change
change_df = pd.DataFrame({
    'Biomolecule': speciesNames,
    'Baseline Steady State': ss,
    'Perturbed Steady State': ss_pertubation,
    'Change in Activity': change_in_activity,
    'Magnitude': np.abs(change_in_activity)  # absolute magnitude for ranking
})

# Sort by absolute change to find the most impacted biomolecules
change_df_sorted = change_df.sort_values(by='Magnitude', ascending=False)

# Display or save the results for one knockdown
print("Most Influential Biomolecules on CCM Phenotypes (sorted by change in activity):")
print(change_df_sorted[['Biomolecule', 'Change in Activity', 'Magnitude']])

# Save the sorted results to a CSV file for further analysis or reporting
#change_df_sorted.to_csv("influential_biomolecules.csv", index=False)

# change_in_activity_ang_PI3k = ss_pertubation[36] - ss[36]
# print(change_in_activity_ang_PI3k)

print("all drugs now")
import matplotlib.pyplot as plt

angiogenesis = 36
change_in_angiogenesis = {}

#loop through all genes except for outputs
for i in range(len(speciesNames) - 3):
    # Copy ymax and set the current biomolecule's ymax to 0 (simulate knockdown)
    #print(speciesNames[i])
    ymax_knockdown = ymax.copy()
    ymax_knockdown[i] = 0

    # Run the simulation with the knockdown
    sol2_dict = {}
    t2_list = []
    y2_list = []
    sol2 = scipy.integrate.solve_ivp(file.ODEfunc,tspan,y0,args=(tau,ymax_knockdown,w,n,EC50,))
    sol2_dict[1] = sol2
    t2_list.append(sol2.t)
    y2_list.append(sol2.y)

    t2 = np.hstack(t2_list)
    y2 = np.hstack(y2_list)

    y2_plot = y2.T

    ss_pertubation = y2_plot[-1,:]  # steady-state values at the end of the simulation

    # Calculate the change in activity for the tracked biomolecule due to the knockdown
    change_in_activity = ss_pertubation[angiogenesis] - ss[angiogenesis]

    # Store the result in the dictionary with the gene name as the key
    change_in_angiogenesis[speciesNames[i]] = change_in_activity

sorted_change = dict(sorted(change_in_angiogenesis.items(), key=lambda item: item[1]))
sorted_change_df = pd.DataFrame(list(sorted_change.items()), columns=['Angiogenesis', 'Change in Activity of Angiogenesis'])
print(sorted_change_df)

plt.figure(figsize=(10, 6))
plt.bar(sorted_change_df['Angiogenesis'], sorted_change_df['Change in Activity of Angiogenesis'])
plt.xticks(rotation=90)
plt.xlabel('Angiogenesis')
plt.ylabel('Change in Activity of Angiogenesis')
plt.title(f'Change in Activity of {speciesNames[angiogenesis]} across Knockdowns')
plt.tight_layout()
plt.show()

# Write to a text file in the same folder
# Convert the array to a space-delimited string
str_change_in_activity = ' '.join(map(str, change_in_activity))
with open("change_in_activity.txt", "w") as file:
    file.write(str_change_in_activity)


#UNIFORM COST SEARCH

#first create the graph

df = pd.read_excel("reactions.xlsx")

# Inspect the first few rows to ensure it's loaded correctly
print(df.head())