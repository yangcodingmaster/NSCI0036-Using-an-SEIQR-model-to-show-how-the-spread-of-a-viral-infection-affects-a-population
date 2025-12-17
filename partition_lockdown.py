import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

from SEIQR_Simulation import update_grid 

# define state constants
Empty = 0
S = 1  # Susceptible
E = 2  # Exposed
I = 3  # Infected
Q = 4  # Quarantined
R = 5  # Recovered

# calculate infection rates for each block and check if lockdown is needed
def update_grid_with_lockdown(grid, time_grid, params, block_indices, 
                             current_step, lockdown_status, threshold=0.1):
    
    # calculate infection rates for each block
    block_infection_rates = []
    for b_idx, indices in enumerate(block_indices):
        # calculate population of the block
        b_population = sum(1 for r, c in indices if grid[r, c] != Empty)
        if b_population == 0:
            block_infection_rates.append(0)
            continue
            
        # calculate infection rate
        infected_count = sum(1 for r, c in indices if grid[r, c] == I)
        infection_rate = infected_count / b_population
        block_infection_rates.append(infection_rate)
        
        # check if lockdown is needed
        if infection_rate >= threshold and not lockdown_status[b_idx]['active']:
            # activate lockdown
            lockdown_status[b_idx]['active'] = True
            lockdown_status[b_idx]['start_time'] = current_step
            print(f"day {current_step}: Block {b_idx} has infection rate {infection_rate:.4f},lockdown activated!")
    
    # store original d value
    original_d = params['d']
    
    # create block parameters dictionary, set d value for each block
    block_params = {}
    for b_idx in range(len(block_indices)):
        block_params[b_idx] = params.copy()
        if lockdown_status[b_idx]['active']:
            block_params[b_idx]['d'] = 0  # set d=0 for locked down blocks
    
    # create cell to block mapping
    L = grid.shape[0]
    cell_to_block = np.full((L, L), -1, dtype=int)  # -1 means not belong to any block
    
    for b_idx, indices in enumerate(block_indices):
        for r, c in indices:
            cell_to_block[r, c] = b_idx
    
    # use the original update_grid method to update the entire grid, the d parameter within the block remains unchanged
    new_grid, new_time_grid = update_grid(grid, time_grid, params)
    
    # apply the update to each locked down block by using d=0
    for b_idx, status in lockdown_status.items():
        if status['active']:
            # create temporary grid and time grid for the block
            temp_grid = grid.copy()
            temp_time_grid = time_grid.copy()
            
            # update the block with d=0
            temp_params = params.copy()
            temp_params['d'] = 0
            
            updated_temp_grid, updated_temp_time_grid = update_grid(temp_grid, temp_time_grid, temp_params)
            
            # merge the update result of the block to the new grid
            for r, c in block_indices[b_idx]:
                new_grid[r, c] = updated_temp_grid[r, c]
                new_time_grid[r, c] = updated_temp_time_grid[r, c]
    
    return new_grid, new_time_grid, params, lockdown_status

# run SEIQR model simulation with lockdown when infection rate exceeds threshold
def run_simulation_with_lockdown(grid, block_indices, params, max_steps=30, threshold=0.1):
    
    L = grid.shape[0]
    time_grid = np.zeros((L, L), dtype=int)
    
    # initialize lockdown status 
    lockdown_status = {}
    for b_idx in range(len(block_indices)):
        lockdown_status[b_idx] = {'active': False, 'start_time': -1}
    
    # store lockdown events
    lockdown_events = []
    
    # store overall data and data for each block
    overall_data = {
        'time': list(range(max_steps)),
        S: [], E: [], I: [], Q: [], R: []
    }
    
    block_data = []
    for _ in range(len(block_indices)):
        block_data.append({
            'time': list(range(max_steps)),
            S: [], E: [], I: [], Q: [], R: []
        })
    
    # store history of grid states
    history = [grid.copy()]
    history_steps = [0]
    
    # main loop
    for step in range(max_steps):
        # update grid and check if lockdown is needed
        grid, time_grid, params, lockdown_status = update_grid_with_lockdown(
            grid, time_grid, params, block_indices, step, lockdown_status, threshold
        )
        
        # collect lockdown events
        for b_idx, status in lockdown_status.items():
            if status['active'] and status['start_time'] == step:
                lockdown_events.append({
                    'block': b_idx,
                    'time': step
                })
        
        # store history every 10 steps
        if (step + 1) % 10 == 0:
            history.append(grid.copy())
            history_steps.append(step + 1)
        
        # calculate overall state 
        total_population = np.sum(grid != Empty)
        if total_population > 0:
            for state in [S, E, I, Q, R]:
                count = np.sum(grid == state)
                overall_data[state].append(count / total_population)
        else:
            for state in [S, E, I, Q, R]:
                overall_data[state].append(0)
        
        # calculate state data for each block
        for b_idx, indices in enumerate(block_indices):
            # Extract the population count of the region
            b_population = sum(1 for r, c in indices if grid[r, c] != Empty)
            
            if b_population > 0:
                for state in [S, E, I, Q, R]:
                    count = sum(1 for r, c in indices if grid[r, c] == state)
                    block_data[b_idx][state].append(count / b_population)
            else:
                for state in [S, E, I, Q, R]:
                    block_data[b_idx][state].append(0)
    
    return overall_data, block_data, lockdown_events, history, history_steps

# plot infection curves for each block with lockdown events
def plot_infection_curves_with_lockdown(block_data, lockdown_events, lockdown_threshold, 
                                      colors=None, save_path=None, show_plot=True):
  
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.figure(figsize=(12, 8))
    
    # set default colors
    if colors is None:
        colors = ['blue', 'green', 'red', 'purple']
    
    # make sure we have enough colors
    if len(colors) < len(block_data):
        # generate default colors if colours are not enough
        from matplotlib.cm import get_cmap
        cmap = get_cmap('tab10')
        colors = [cmap(i) for i in range(len(block_data))]
    
    for b_idx, data in enumerate(block_data):
        time_steps = data['time']
        infected = data[I]
        
        # plot infection curve
        plt.plot(time_steps, infected, color=colors[b_idx], linewidth=2,
                label=f'Block {b_idx + 1}')
        
        # mark the maximum infection rate point
        max_idx = np.argmax(infected)
        max_val = infected[max_idx]
        plt.scatter(time_steps[max_idx], max_val, color=colors[b_idx], s=80, zorder=5, edgecolor='black')
        plt.annotate(f'{max_val:.3f}', xy=(time_steps[max_idx], max_val), xytext=(5, 5),
                    textcoords='offset points', fontsize=8)
        
        # mark the lockdown time
        for event in lockdown_events:
            if event['block'] == b_idx:
                lockdown_time = event['time']
                lockdown_rate = infected[lockdown_time]
                plt.axvline(x=lockdown_time, color=colors[b_idx], linestyle='--', alpha=0.5)
                plt.scatter(lockdown_time, lockdown_rate, color=colors[b_idx], marker='X', s=100, zorder=5)
                plt.annotate(f'Lockdown', xy=(lockdown_time, lockdown_rate), 
                            xytext=(10, -20), textcoords='offset points',
                            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'),
                            fontsize=8)
    
    # draw threshold line
    plt.axhline(y=lockdown_threshold, color='black', linestyle=':', alpha=0.7,
               label=f'Lockdown Threshold ({lockdown_threshold*100}%)')
    
    # set plot properties
    plt.xlabel('day', fontsize=15)
    plt.ylabel('Infected Population Fraction', fontsize=15)
    #plt.title('Infection Curves by Block with Lockdown', fontsize=15)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # save the plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # show the plot
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return plt

# compare infection curves with and without lockdown for each block
def plot_lockdown_comparison(block_data, block_data_no_lockdown, lockdown_events, 
                           lockdown_threshold, block_indices=None, block_name="Block",
                           layout=None, figsize=None, save_path=None, show_plot=True):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    
    # confirmed the regions to be plotted
    if block_indices is None:
        # if not specified, plot all regions
        block_indices = list(range(len(block_data)))
    
    # number of blocks
    block_count = len(block_indices)
    
    # calculate subplot layout
    if layout is None:
        if block_count <= 1:
            rows, cols = 1, 1
        elif block_count <= 2:
            rows, cols = 1, 2
        elif block_count <= 4:
            rows, cols = 2, 2
        else:
            # calculate suitable layout for more blocks
            cols = math.ceil(math.sqrt(block_count))
            rows = math.ceil(block_count / cols)
    else:
        rows, cols = layout
    
    # set figure size
    if figsize is None:
        figsize = (5*cols, 5*rows)
    
    # create figure
    plt.figure(figsize=figsize)
    
    # plot each specified block
    for i, b_idx in enumerate(block_indices):
        plt.subplot(rows, cols, i+1)
        
        # with lockdown
        time_steps = block_data[b_idx]['time']
        infected_with_lockdown = block_data[b_idx][I]
        plt.plot(time_steps, infected_with_lockdown, 'r-', label='With Lockdown', linewidth=2)
        
        # without lockdown
        infected_no_lockdown = block_data_no_lockdown[b_idx][I]
        plt.plot(time_steps, infected_no_lockdown, 'b--', label='Without Lockdown', linewidth=2)
        
        # mark the lockdown time for the block
        for event in lockdown_events:
            if event['block'] == b_idx:
                lockdown_time = event['time']
                plt.axvline(x=lockdown_time, color='red', linestyle='--', alpha=0.5)
                plt.text(lockdown_time, 0.01, 'Lockdown', 
                        rotation=90, verticalalignment='bottom', fontsize=8)
        
        # draw threshold line
        plt.axhline(y=lockdown_threshold, color='black', linestyle=':', alpha=0.7)
        
        # set plot properties
        plt.xlabel('day', fontsize=15)
        plt.ylabel('Infected Population Fraction', fontsize=15)
        #plt.title(f'{block_name} {b_idx + 1} Comparison', fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # ensure y-axis starts from 0
        plt.ylim(bottom=0)
        
        # if there is data, find the maximum value and set the y-axis limit
        if len(infected_with_lockdown) > 0 and len(infected_no_lockdown) > 0:
            max_val = max(max(infected_with_lockdown), max(infected_no_lockdown))
            # add a little margin
            plt.ylim(top=min(1.0, max_val * 1.1))
    
    plt.tight_layout()
    
    # save the plot
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    # show the plot
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return plt
