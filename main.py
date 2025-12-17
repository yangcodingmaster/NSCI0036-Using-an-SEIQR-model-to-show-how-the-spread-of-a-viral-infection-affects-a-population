from SEIQR_Simulation import initialised_grid, update_grid, calculate_state_proportions, run_simulation, parameter_scan
from SEIQR_Plot import plot_parameter_impact_on_infection, plot_parameter_impact_on_quarantine, plot_max_infection_vs_parameter, plot_all_states, plot_all_states_for_parameter, visualize_spatial_spread
from partition_simulation import initialised_partitioned_grid, update_grid, run_partitioned_simulation, visualize_partitioned_grid, plot_block_infection_curves, visualize_simulation_history
from partition_lockdown import update_grid_with_lockdown, run_simulation_with_lockdown, plot_infection_curves_with_lockdown, plot_lockdown_comparison

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import time

# Define state constants
Empty = 0
S = 1  # Susceptible
E = 2  # Exposed
I = 3  # Infected
Q = 4  # Quarantined
R = 5  # Recovered

# set base parameters (default)
base_params = {
    'pe': 0.5,    # probability of exposure
    'pi': 0.5,    # probability of infection
    'pq': 0.1,    # probability of quarantine
    'pr': 0.12,   # prabability of recovery
    'ti': 7,      # incubation
    'tq': 7,      # time taken for detection
    'tr': 21,     # time taken for recovery
    'd': 2        # Moore neighbourhood
}

# pe
pe_values = [0.1, 0.3, 0.5, 0.7, 0.9]
pe_results, pe_metadata = parameter_scan('pe', pe_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_parameter_impact_on_infection(pe_results, pe_metadata, 
                                 save_path='pe_impact_on_infection.png')

# pi
pi_values = [0.01, 0.05, 0.1, 0.5, 1]
pi_results, pi_metadata = parameter_scan('pi', pi_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_parameter_impact_on_infection(pi_results, pi_metadata, 
                                 save_path='pi_impact_on_infection.png')

# pq
pq_values = [0.01, 0.05, 0.1, 0.5, 1]
pq_results, pq_metadata = parameter_scan('pq', pq_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_parameter_impact_on_infection(pq_results, pq_metadata, 
                                 save_path='pq_impact_on_infection.png')

pq_values = [0.01, 0.05, 0.1, 0.5, 1]
pq_results, pq_metadata = parameter_scan('pq', pq_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_parameter_impact_on_quarantine(pq_results, pq_metadata, 
                                 save_path='pq_impact_on_quarantine.png')

# pr
pr_values = [0.001, 0.01, 0.1, 0.5, 1]
pr_results, pr_metadata = parameter_scan('pr', pr_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_parameter_impact_on_infection(pr_results, pr_metadata, 
                                 save_path='pr_impact_on_infection.png', xlim=(20, 80))

# relation between max infection fraction and pi
pi_values = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1]
pi_results, pi_metadata = parameter_scan('pi', pi_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_max_infection_vs_parameter(pi_metadata, title='i_max agaist time of different p_i values', save_path='pi_imax.png')

# relation betweem max infection fraction and pq
pq_values = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1]
pq_results, pq_metadata = parameter_scan('pq', pq_values, base_params, 
                                        grid_size=100, 
                                        population_density=0.75,
                                        init_exposure=50, 
                                        init_infection=5,
                                        max_steps=100)
plot_max_infection_vs_parameter(pq_metadata, title='i_max agaist time of different p_q values', save_path='pq_imax.png')

# ti effect of incubation
ti_values = [1, 7, 14, 21, 28]
ti_results, ti_metadata = parameter_scan('ti', ti_values, base_params,
                                         max_steps= 100)
plot_parameter_impact_on_infection(ti_results, ti_metadata, 
                                 save_path='ti_impact_on_infection.png')

# tr effect of recoverey
tr_values = [1, 7, 14, 21, 28]
tr_results, tr_metadata = parameter_scan('tr', tr_values, base_params,
                                         max_steps= 100)
plot_parameter_impact_on_infection(tr_results, tr_metadata, 
                                 save_path='tr_impact_on_infection.png', xlim=(20, 80))
plot_parameter_impact_on_quarantine(tr_results, tr_metadata, 
                                 save_path='tr_impact_on_quarantine.png')

# tq effect of recoverey
tq_values = [1, 7, 14, 21, 28]
tq_results, tq_metadata = parameter_scan('tq', tr_values, base_params,
                                         max_steps= 100)
plot_parameter_impact_on_infection(tq_results, tq_metadata, 
                                 save_path='tq_impact_on_infection.png', xlim=(20,80))
plot_parameter_impact_on_quarantine(tq_results, tq_metadata, 
                                 save_path='tq_impact_on_quarantine.png')

# default parameters impact on all state curve
# run a simulation with default parameters
default_params = {
    'pe': 0.5,    # probability of exposure
    'pi': 0.5,    # probability of infection
    'pq': 0.1,    # probability of quarantine
    'pr': 0.12,   # prabability of recovery
    'ti': 7,      # incubation
    'tq': 7,      # time taken for detection
    'tr': 21,     # time taken for recovery
    'd': 3        # Moore neighbourhood
}

# initialisation
grid_size = 100
grid = initialised_grid(grid_size, population_density=0.75, init_Exposure=50, init_Infection=5)

# run a simulation
default_data = run_simulation(grid, default_params, max_steps=100)

# plot all states
plot_all_states(default_data, params=default_params, save_path='default_simulation.png', xlim=(0, 100))

# spread simulation
d_values = [1, 2, 3]
d_results, d_metadata = parameter_scan('d', d_values, base_params)

for d_value in d_values:
    params= base_params.copy()
    params['d'] = d_value
    visualize_spatial_spread(L=100,
                             params=params,
                             total_steps=60,
                             save_interval=20,
                             figsize=(12, 4), save_path=f'SEIQR_model_simulation{d_value}.png')

# set parameters
L = 100  # grid size (100x100)
population_density = 0.75  # population density
init_Exposure = 50  # initial exposed
init_Infection = 5  # initial infected
infected_block = 0  # initial infected block left top corner

    
# initialisation grid with cross-shaped Empty regions
print("initialising grid with cross-shaped Empty regions...")
start_time = time.time()
grid, block_indices = initialised_partitioned_grid(
    L, population_density, init_Exposure, init_Infection, infected_block
)
print(f"finish initialisation, time taken: {time.time() - start_time:.2f} seconds")


# set lockdown threshold
lockdown_threshold = 0.03  # lockdown when infection rate reaches 3%
lockdown_thresholds = [0.01, 0.05, 0.10, 0.15, 0.20]

for threshold in lockdown_thresholds:
    # run simulation with lockdown
    print(f"run simulation with infection threshold {threshold*100}%...")
    start_time = time.time()
    max_steps = 30
    overall_data, block_data, lockdown_events, history, history_steps = run_simulation_with_lockdown(
        grid, block_indices, params, max_steps, threshold
    )
    print(f"finish simulation, time taken: {time.time() - start_time:.2f} seconds")

    # output lockdown events
    if lockdown_events:
        print("\nrecord the lockdown event:")
        for event in lockdown_events:
            print(f"day {event['time']}: block {event['block'] + 1} is locked down")
    else:
        print("\nno lockdown event recorded")

    # visualise simulation history
    print("visualise simulation history...")
    visualize_simulation_history(history, history_steps, block_indices, save_path="block_simulation.png")



    plot_infection_curves_with_lockdown(
        block_data, 
        lockdown_events, 
        threshold,
        save_path="infection_curves_with_lockdown.png"
    )

    # compare the results with and without lockdown
    print("comparing with and without lockdown...")

    # running simulation without lockdown
    print("running simulation without lockdown...")
    start_time = time.time()
    grid_no_lockdown, _ = initialised_partitioned_grid(
        L, population_density, init_Exposure, init_Infection, infected_block
    )
    overall_data_no_lockdown, block_data_no_lockdown, _, _ = run_partitioned_simulation(
        grid_no_lockdown, block_indices, params, max_steps
    )
    print(f"comparison finished, time taken: {time.time() - start_time:.2f}seconds")


    plot_lockdown_comparison(
        block_data, 
        block_data_no_lockdown, 
        lockdown_events, 
        threshold,
        block_indices=[0, 3], # compare the block 1 and block 4
        layout=(1, 2),  #one row two columns
        save_path=f"specific_blocks_comparison{threshold}.png"
    )




