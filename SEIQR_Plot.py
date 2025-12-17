import matplotlib.colors as mcolors
import time

from SEIQR_Simulation import initialised_grid, update_grid 

# define state constants
Empty = 0
S = 1  # Susceptible
E = 2  # Exposed
I = 3  # Infected
Q = 4  # Quarantined
R = 5  # Recovered

# plot parameter impact on infection
def plot_parameter_impact_on_infection(results, metadata, save_path=None, show_max_points=True, xlim=None):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    param_name = metadata['param_name']
    param_values = list(results.keys())
    
    plt.figure(figsize=(12, 8))
    
    # use different colors and markers
    colors = plt.cm.get_cmap('viridis')(np.linspace(0, 1, len(param_values)))
    #markers = ['o', 's', '^', 'D', 'x', 'p', '*', '+']
    
    for i, value in enumerate(param_values):
        data = results[value]
        time = data['time']
        infected = data[I]
        
        color = colors[i]
        #marker = markers[i % len(markers)]
        
        plt.plot(time, infected, color=color, 
                 marker=None, markevery=max(1, len(time)//10),
                 label=f'{param_name} = {value}', linewidth=2, markersize=8)
        
        # mark the max infection point
        if show_max_points:
            max_idx = np.argmax(infected)
            max_val = infected[max_idx]
            plt.scatter(time[max_idx], max_val, color=color, 
                       s=100, zorder=5, edgecolor='black')
    
    # set plot properties
    plt.xlabel('Day', fontsize=15)
    plt.ylabel('infected population fraction', fontsize=15)
    #plt.title(f'effect of parameter {param_name}', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    # set y axis range
    plt.ylim(bottom=0)
    
     # set x axis limit if have
    if xlim is not None:
        plt.xlim(xlim)
    
    # add grid lines
    #plt.grid(True, linestyle='--', alpha=0.7)
    
    # plot comments
    base_params_text = ', '.join([f"{k}={v}" for k, v in metadata['base_params'].items() if k != param_name])
    #plt.figtext(0.01, 0.01, f"fixed parameters: {base_params_text}", 
                #ha="left", fontsize=10, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    plt.tight_layout()
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()
 
# plot parameter impact on quarantine    
def plot_parameter_impact_on_quarantine(results, metadata, save_path=None, show_max_points=True, xlim=None):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    param_name = metadata['param_name']
    param_values = list(results.keys())
    
    plt.figure(figsize=(12, 8))
    
    # use different colors and markers
    colors = plt.cm.get_cmap('viridis')(np.linspace(0, 1, len(param_values)))
    #markers = ['o', 's', '^', 'D', 'x', 'p', '*', '+']
    
    for i, value in enumerate(param_values):
        data = results[value]
        time = data['time']
        qurantine = data[Q]
        
        color = colors[i]
        #marker = markers[i % len(markers)]
        
        plt.plot(time, qurantine, color=color, 
                 marker=None, markevery=max(1, len(time)//10),
                 label=f'{param_name} = {value}', linewidth=2, markersize=8)
        
        # mark the max infection point
        if show_max_points:
            max_idx = np.argmax(qurantine)
            max_val = qurantine[max_idx]
            plt.scatter(time[max_idx], max_val, color=color, 
                       s=100, zorder=5, edgecolor='black')
    
    # set plot properties
    plt.xlabel('Day', fontsize=15)
    plt.ylabel('qurantine population fraction', fontsize=15)
    #plt.title(f'effect of parameter {param_name}', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    # set y axis range
    plt.ylim(bottom=0)
    
     # set x axis limit if have
    if xlim is not None:
        plt.xlim(xlim)
    
    # add grid lines
    #plt.grid(True, linestyle='--', alpha=0.7)
    
    # plot comments
    base_params_text = ', '.join([f"{k}={v}" for k, v in metadata['base_params'].items() if k != param_name])
    #plt.figtext(0.01, 0.01, f"fixed parameters: {base_params_text}", 
                #ha="left", fontsize=10, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    plt.tight_layout()
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# plot the maximum infection rate for a parameter value
def plot_max_infection_vs_parameter(metadata, title=None, save_path=None):
   
    import matplotlib.pyplot as plt
    import numpy as np
    
    param_name = metadata['param_name']
    param_values = sorted(list(metadata['max_infected'].keys()))
    max_infected_values = [metadata['max_infected'][val] for val in param_values]
    
    plt.figure(figsize=(10, 6))
    
    # relation between max infection fraction and a paramter
    plt.plot(param_values, max_infected_values, 'o-', color='crimson', 
             linewidth=2, markersize=10)
    
    # best fit the trend
    if len(param_values) > 2:  # at least 3 points required to fit
        try:
            z = np.polyfit(param_values, max_infected_values, 2)
            p = np.poly1d(z)
            x_smooth = np.linspace(min(param_values), max(param_values), 100)
            plt.plot(x_smooth, p(x_smooth), '--', color='gray', alpha=0.7)
        except:
            pass  # pass if fitting failed
    
    # mark the data point
    for x, y in zip(param_values, max_infected_values):
        plt.annotate(f'{y:.3f}', 
                    xy=(x, y), 
                    xytext=(0, 10),
                    textcoords='offset points',
                    ha='center', 
                    fontsize=9)
    
    # set plot properties
    plt.xlabel(f'value of parameter {param_name} ', fontsize=15)
    plt.ylabel('infected population fraction', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    #if title:
        #plt.title(title, fontsize=16)
    #else:
        #plt.title(f'Effect of paramter {param_name} to max infection rate', fontsize=16)
    
    #plt.grid(True, linestyle='--', alpha=0.7)
    
    # fixed paramters comment
    base_params_text = ', '.join([f"{k}={v}" for k, v in metadata['base_params'].items() if k != param_name])
    #plt.figtext(0.01, 0.01, f"fixed: {base_params_text}", 
                #ha="left", fontsize=10, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    plt.tight_layout()
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# plot all the states that infection rate against time
def plot_all_states(data, params=None, title="SEIQR Model Default Simulation", save_path=None, xlim=None):
   
    import matplotlib.pyplot as plt
    
    time = data['time']
    
    plt.figure(figsize=(12, 7))
    
    # plot state curves
    plt.plot(time, data[S], 'b-', label='susceptible (S)', linewidth=2)
    plt.plot(time, data[E], 'g-', label='exposed (E)', linewidth=2)
    plt.plot(time, data[I], 'r-', label='infected (I)', linewidth=2)
    plt.plot(time, data[Q], 'c-', label='quarantined (Q)', linewidth=2)
    plt.plot(time, data[R], 'y-', label='recovered (R)', linewidth=2)
    
    # add grey dash grid
    plt.grid(True, linestyle='--', alpha=0.7, color='gray')
    
    # set plot properties
    plt.xlabel('days', fontsize=15)
    plt.ylabel('population fraction', fontsize=15)
    #plt.title(title, fontsize=14)
    plt.legend(loc='center right', fontsize=10)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    # set x axis limit if have
    if xlim is not None:
        plt.xlim(xlim)
    
    # find the corresponding time of max infection fraction
    max_infected_idx = data[I].index(max(data[I]))
    max_infected = data[I][max_infected_idx]
    max_infected_time = time[max_infected_idx]
    
    # mark the  max infection point
    plt.scatter(max_infected_time, max_infected, color='darkred', s=100, zorder=5)
    plt.annotate(f'max infection rate: {max_infected:.3f}', 
                xy=(max_infected_time, max_infected),
                xytext=(max_infected_time + 5, max_infected),
                fontsize=10,
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5))
    
    # show parameters at bottom left if provided
    if params:
        params_text = ', '.join([f"{k}={v}" for k, v in params.items()])
        #plt.figtext(0.01, 0.01, f"Parameters: {params_text}", 
                   #ha="left", fontsize=9, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    plt.tight_layout()
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# plot all the statees for a parameter value
def plot_all_states_for_parameter(results, param_value, metadata, save_path=None, xlim=None):
  
    import matplotlib.pyplot as plt
    
    param_name = metadata['param_name']
    
    # make sure the parameter value is in the result
    if param_value not in results:
        available_values = list(results.keys())
        param_value = min(available_values, key=lambda x: abs(x - param_value))
        print(f"warning: parameter {param_value} is not in the result, using the closest value: {param_value}")
    
    data = results[param_value]
    time = data['time']
    
    plt.figure(figsize=(12, 7))
    
    # plot state curves
    plt.plot(time, data[S], 'b-', label='susceptible (S)', linewidth=2)
    plt.plot(time, data[E], 'g-', label='exposed (E)', linewidth=2)
    plt.plot(time, data[I], 'r-', label='infected (I)', linewidth=2)
    plt.plot(time, data[Q], 'c-', label='quarantined (Q)', linewidth=2)
    plt.plot(time, data[R], 'y-', label='recovered (R)', linewidth=2)
    
    # dash grey grid
    plt.grid(True, linestyle='--', alpha=0.7, color='gray')
    
    # set plot properties
    plt.xlabel('days', fontsize=15)
    plt.ylabel('population fraction', fontsize=15)
    #plt.title(f'SEIQR model simulation curve ({param_name} = {param_value})', fontsize=14)
    plt.legend(loc='center right', fontsize=12)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    # set x axis limit if have
    if xlim is not None:
        plt.xlim(xlim)
    
    # find the corresponding time of max infection fraction
    max_infected_idx = data[I].index(max(data[I]))
    max_infected = data[I][max_infected_idx]
    max_infected_time = time[max_infected_idx]
    
    # mark the max infection 
    plt.scatter(max_infected_time, max_infected, color='darkred', s=100, zorder=5)
    plt.annotate(f'max infection rate: {max_infected:.3f}', 
                xy=(max_infected_time, max_infected),
                xytext=(max_infected_time + 5, max_infected),
                fontsize=10,
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5))
    
    plt.tight_layout()
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# Visualisation of spatial spread of SEIQR model
def visualize_spatial_spread(L=200, population_density=0.5, init_exposure=50, init_infection=5, 
                           params=None, total_steps=100, save_interval=20, figsize=(15, 6), save_path = None):
   
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches
    
    
    # default paramters
    if params is None:
        params = {
            'pe': 0.5, 'pi': 0.5, 'pq': 0.1, 'pr': 0.12,
            'ti': 7, 'tq': 7, 'tr': 21, 'd': 2
        }
    
    # initailisation
    current_grid = initialised_grid(L, population_density, init_exposure, init_infection)
    current_time = np.zeros((L, L), dtype=int)
    
    # store state history 
    history = [current_grid.copy()]
    
    # run simulation
    for step in range(total_steps):
        current_grid, current_time = update_grid(current_grid, current_time, params)
        if (step + 1) % save_interval == 0:
            history.append(current_grid.copy())
    
    # visualisation
    n_images = len(history)
    fig, axes = plt.subplots(1, n_images, figsize=figsize)
    plt.subplots_adjust(top=0.85, bottom=0.15, hspace=0.3, wspace=0.1)
    
    # creating colour mapping
    cmap = mcolors.ListedColormap(['grey', 'purple', 'orange', 'red', 'blue', 'green'])
    state_labels = ['Empty', 'Susceptible', 'Exposed', 'Infected', 'Quarantined', 'Recovered']
    
    # draw state at each timepstep
    for idx, ax in enumerate(axes.flatten()):
        if idx >= len(history): 
            break
        im = ax.imshow(history[idx], cmap=cmap, vmin=0, vmax=5, interpolation='nearest')
        ax.set_title(f"day {idx*save_interval}", fontsize=15)
        ax.axis('off')
    
    # legend
    patches = [mpatches.Patch(color=cmap(i), label=label) 
              for i, label in enumerate(state_labels)]
    
    fig.legend(
        handles=patches,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.1),  
        ncol=6,
        fontsize=12,
        frameon=False
    )
    
    # title
    param_str = f"pe={params['pe']}, pi={params['pi']}, pq={params['pq']}, pr={params['pr']}, d={params['d']}"
    #plt.suptitle(f"SEIQR Model Spatial Spread Simulation", y=0.85, fontsize=12)
    
    #plt.figtext(0.5, 0.24, f"parameters: {param_str}", 
                #ha="center", fontsize=8, bbox={"facecolor":"orange", "alpha":0.2, "pad":4})
    
    # save plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()
    
    return history

