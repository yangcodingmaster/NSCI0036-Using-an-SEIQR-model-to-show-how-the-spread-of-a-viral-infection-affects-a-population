import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# define state constants
Empty = 0
S = 1  # Susceptible
E = 2  # Exposed
I = 3  # Infected
Q = 4  # Quarantined
R = 5  # Recovered

# initialise a grid and use a cross-shaped Empty region with width 2 to partition the grid into four blocks
def initialised_partitioned_grid(L, population_density, init_Exposure, init_Infection, infected_block=0):
    
    # make sure L is even
    if L % 2 != 0:
        L += 1
        print(f"grid size length with even number: {L}")
    
    # create an empty grid
    grid = np.zeros((L, L), dtype=int)
    
    # define the width of the cross
    cross_width = 2  # 2 cells wide
    
    # calculate the center of the grid
    center = L // 2
    
    # create a cross-shaped empty region
    # horizontal line
    grid[center-cross_width//2:center+cross_width//2, :] = Empty
    # vertical line
    grid[:, center-cross_width//2:center+cross_width//2] = Empty
    
    # define the boundaries of each block
    # Note: boundaries should consider the cross region to avoid overlap
    block_boundaries = [
        (0, center-cross_width//2, 0, center-cross_width//2),                # left up: (min_row, max_row, min_col, max_col)
        (0, center-cross_width//2, center+cross_width//2, L),                # right up
        (center+cross_width//2, L, 0, center-cross_width//2),                # left down
        (center+cross_width//2, L, center+cross_width//2, L)                 # right down
    ]
    
    # create a list to store the indices of each block
    block_indices = []
    for min_row, max_row, min_col, max_col in block_boundaries:
        indices = []
        for r in range(min_row, max_row):
            for c in range(min_col, max_col):
                indices.append((r, c))
        block_indices.append(indices)
    
    # randomly assign population to each block
    for b_idx, indices in enumerate(block_indices):
        # calculate the total number of cells
        total_cells = len(indices)
        population = int(population_density * total_cells)
        
        # randomly shuffle the indices to allocate population
        np.random.shuffle(indices)
        occupied_cells = indices[:population]
        
        # assign the S state to the occupied cells
        for r, c in occupied_cells:
            grid[r, c] = S
    
    # assign initial exposure and initial infection at pointed block
    infected_indices = block_indices[infected_block]
    S_positions = [(r, c) for r, c in infected_indices if grid[r, c] == S]
    np.random.shuffle(S_positions)
    
    # assign initial exposure
    for i in range(min(init_Exposure, len(S_positions))):
        r, c = S_positions[i]
        grid[r, c] = E
    
    # update available S positions
    S_positions = [(r, c) for r, c in infected_indices if grid[r, c] == S]
    np.random.shuffle(S_positions)
    
    # assign initial infection
    for i in range(min(init_Infection, len(S_positions))):
        r, c = S_positions[i]
        grid[r, c] = I
    
    return grid, block_indices

# update grid in every timestep
def update_grid(grid, time_grid, params):
    L = grid.shape[0]
    
    d = params['d']
    neighbor_counts = [(dx, dy) for dx in range(-d, d+1)
                       for dy in range(-d, d+1)
                       if not (dx == 0 and dy == 0)
    ]
    
    # store states and state durations
    new_grid = grid.copy()
    new_time_grid = time_grid.copy()

    # loop every cell, including boundary
    for i in range(L):
        for j in range(L):
            current_state = grid[i, j]
            # do not update if the cell is recovered or empty
            if current_state == Empty or current_state == R:
                continue
            
            # current cell state duration
            duration = time_grid[i, j]
            
            # update the S to E without time delay
            if current_state == S:
                # count number of E and I in neighbor
                infected_neighbors = 0
                for dx, dy in neighbor_counts:
                    ni, nj = i + dx, j + dy
                    if 0 <= ni < L and 0 <= nj < L:
                        if grid[ni, nj] in [E, I]:
                            infected_neighbors += 1
                # state conversion descision by p and number of E and I
                if infected_neighbors > 0 and np.random.rand() < params["pe"]:
                    new_grid[i, j] = E
                    new_time_grid[i, j] = 0  # reset duration
                else:
                    new_time_grid[i, j] += 1
            
            # update state E
            elif current_state == E:
                # state converts to I if t > ti
                if duration >= params['ti']:
                    if np.random.rand() < params['pi']:
                        new_grid[i, j] = I
                        new_time_grid[i, j] = 0  # reset duration
                    else:
                        # stay E
                        new_time_grid[i, j] += 1
                else:
                    new_time_grid[i, j] += 1
            
            # update I to Q or R
            elif current_state == I:
                # converts to Q based on pq and time
                if duration >= params['tq'] and np.random.rand() < params['pq']:
                    new_grid[i, j] = Q
                    new_time_grid[i, j] = 0
                # or attempt to recover
                elif duration >= params['tr'] and np.random.rand() < params['pr']:
                    new_grid[i, j] = R
                    new_time_grid[i, j] = 0
                else:
                    new_time_grid[i, j] += 1
            
            # update Q to R
            elif current_state == Q:
                if duration >= params['tr']:
                    if np.random.rand() < params['pr']:
                        new_grid[i, j] = R
                        new_time_grid[i, j] = 0
                    else:
                        new_time_grid[i, j] += 1  # stay Q
                else:
                    new_time_grid[i, j] += 1

    return new_grid, new_time_grid

# run simulation with partitioned grid with time steps and collect data for each block
def run_partitioned_simulation(grid, block_indices, params, max_steps=100):
    
    L = grid.shape[0]
    time_grid = np.zeros((L, L), dtype=int)
    
    # store overall data and block data
    overall_data = {
        'time': list(range(max_steps)),
        S: [], E: [], I: [], Q: [], R: []
    }
    
    block_data = []
    for _ in range(4):
        block_data.append({
            'time': list(range(max_steps)),
            S: [], E: [], I: [], Q: [], R: []
        })
    
    # history data for visualization
    history = [grid.copy()]
    history_steps = [0]
    
    # main siimulation loop
    for step in range(max_steps):
        # update grid
        grid, time_grid = update_grid(grid, time_grid, params)
        
        # store state every 10 steps
        if (step + 1) % 10 == 0:
            history.append(grid.copy())
            history_steps.append(step + 1)
        
        # calculate overall data
        total_population = np.sum(grid != Empty)
        if total_population > 0:
            for state in [S, E, I, Q, R]:
                count = np.sum(grid == state)
                overall_data[state].append(count / total_population)
        else:
            for state in [S, E, I, Q, R]:
                overall_data[state].append(0)
        
        # calculate block data
        for b_idx, indices in enumerate(block_indices):
            # calculate the population in the block
            b_population = sum(1 for r, c in indices if grid[r, c] != Empty)
            
            if b_population > 0:
                for state in [S, E, I, Q, R]:
                    count = sum(1 for r, c in indices if grid[r, c] == state)
                    block_data[b_idx][state].append(count / b_population)
            else:
                for state in [S, E, I, Q, R]:
                    block_data[b_idx][state].append(0)
    
    return overall_data, block_data, history, history_steps

# visualize partitioned grid
def visualize_partitioned_grid(grid, block_indices, title=""):

    L = grid.shape[0]
    
    plt.figure(figsize=(10, 10))
    
    # create a color map
    cmap = mcolors.ListedColormap(['grey', 'purple', 'orange', 'red', 'blue', 'green'])
    state_labels = ['Empty', 'Susceptible', 'Exposed', 'Infected', 'Quarantined', 'Recovered']
    
    # plot the grid
    plt.imshow(grid, cmap=cmap, vmin=0, vmax=5, interpolation='nearest')
    
    # add legend
    patches = [mpatches.Patch(color=cmap(i), label=label) 
              for i, label in enumerate(state_labels)]
    plt.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, 1.15), 
              ncol=6, fontsize=10)
    
    # add block labels
    block_labels = ['block 1', 'block 2', 'block 3', 'block 4']
    for b_idx, indices in enumerate(block_indices):
        if indices:  # make sure the block is not empty
            # calculate the center of the block
            center_r = sum(r for r, _ in indices) / len(indices)
            center_c = sum(c for _, c in indices) / len(indices)
            plt.text(center_c, center_r, block_labels[b_idx], 
                    color='black', ha='center', va='center', fontsize=12,
                    bbox=dict(facecolor='white', alpha=0.7))
    
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

# plot infection curves by block
def plot_block_infection_curves(block_data, params=None, title="Infection Curves by Quadrant", save_path=None):

    plt.figure(figsize=(12, 8))
    
    colors = ['blue', 'green', 'red', 'purple']
    
    # plot the infection curves for each block
    for b_idx, data in enumerate(block_data):
        time = data['time']
        infected = data[I]
        
        plt.plot(time, infected, color=colors[b_idx], linewidth=2,
                label=f'Block {b_idx + 1}')
        
        # annotate the peak value
        max_idx = np.argmax(infected)
        max_val = infected[max_idx]
        plt.scatter(time[max_idx], max_val, color=colors[b_idx], s=80, zorder=5, edgecolor='black')
        plt.annotate(f'{max_val:.3f}', xy=(time[max_idx], max_val), xytext=(5, 5),
                    textcoords='offset points', fontsize=8)
    
    # add labels and title
    plt.xlabel('day', fontsize=12)
    plt.ylabel('Infected Population Fraction', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # show parameters at the lower left corner
    if params:
        params_text = ', '.join([f"{k}={v}" for k, v in params.items()])
        plt.figtext(0.01, 0.01, f"Parameters: {params_text}", 
                   ha="left", fontsize=9, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# visualize simulation history
def visualize_simulation_history(history, history_steps, block_indices=None, save_path = None):

    n_frames = len(history)
    
    # create a color map
    cmap = mcolors.ListedColormap(['grey', 'purple', 'orange', 'red', 'blue', 'green'])
    state_labels = ['Empty', 'Susceptible', 'Exposed', 'Infected', 'Quarantined', 'Recovered']
    
    fig, axes = plt.subplots(1, n_frames, figsize=(5*n_frames, 6))
    plt.subplots_adjust(top=0.85, bottom=0.15, hspace=0.3, wspace=0.1)
    
    # plot each frame
    for idx, ax in enumerate(axes.flatten()):
        if idx < n_frames:
            im = ax.imshow(history[idx], cmap=cmap, vmin=0, vmax=5, interpolation='nearest')
            ax.set_title(f"day {history_steps[idx]}", fontsize=20)
            ax.axis('off')
            
            # add block labels if have
            if block_indices:
                for b_idx, indices in enumerate(block_indices):
                    if indices:
                        center_r = sum(r for r, _ in indices) / len(indices)
                        center_c = sum(c for _, c in indices) / len(indices)
                        ax.text(center_c, center_r, f'block {b_idx + 1}', 
                               color='black', ha='center', va='center', fontsize=12,
                               bbox=dict(facecolor='white', alpha=0.9))
    
    # add legend
    patches = [mpatches.Patch(color=cmap(i), label=label) 
              for i, label in enumerate(state_labels)]
    
    fig.legend(
        handles=patches,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.02),
        ncol=6,
        fontsize=15,
        frameon=False
    )
    
    #plt.suptitle("SEIQR Model Spatial Spread Simulation with Partitioned Grid", y=0.95, fontsize=14)
    plt.tight_layout(rect=(0, 0.05, 1, 0.95))
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
    plt.show()
