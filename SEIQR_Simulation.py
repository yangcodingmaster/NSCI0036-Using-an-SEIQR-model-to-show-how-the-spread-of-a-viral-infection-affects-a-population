import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import time

# Define state constants
Empty = 0
S = 1  # Susceptible
E = 2  # Exposed
I = 3  # Infected
Q = 4  # Quarantined
R = 5  # Recovered

# initialise a grid, allocate populations, initial exposure, initial infection
def initialised_grid(L, population_density, init_Exposure, init_Infection):
    grid = np.zeros((L, L), dtype=int)
    
    total_cells = L * L
    population = int(population_density * total_cells)
    
    #allocating population in the city
    all_indices = np.arange(total_cells)
    np.random.shuffle(all_indices)
    occupied_cells = all_indices[:population]
    
    #setting the population to be susceptible, S
    for i in occupied_cells:
        r = i // L
        c = i % L
        grid[r, c] = S
        
    #allocating exposure popluation among suscepitable
    S_positions = np.argwhere(grid == S)
    np.random.shuffle(S_positions)
    
    for i in range(min(init_Exposure, len(S_positions))):
        r, c = S_positions[i]
        grid[r, c] = E
        
    #allocating infected population among susceptible
    S_positions = np.argwhere(grid == S)
    np.random.shuffle(S_positions)
    
    for i in range(min(init_Infection, len(S_positions))):
        r, c = S_positions[i]
        grid[r, c] = I
        
    return grid

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
                        new_time_grid[i, j] += 1  # 处理概率失败的情况
                else:
                    new_time_grid[i, j] += 1

    return new_grid, new_time_grid

# calculate population fraction of each state
def calculate_state_proportions(grid):
    # calculate total population exclude empty cell
    total_population = np.sum(grid != Empty)
    
    if total_population == 0:
        return {S: 0, E: 0, I: 0, Q: 0, R: 0}
    
    # count number of each state
    state_counts = {
        S: np.sum(grid == S),
        E: np.sum(grid == E),
        I: np.sum(grid == I),
        Q: np.sum(grid == Q),
        R: np.sum(grid == R)
    }
    
    # converts to ratio
    state_proportions = {state: count / total_population 
                         for state, count in state_counts.items()}
    
    return state_proportions

# operate one simulation in max_steps days
def run_simulation(grid, params, max_steps=100):
    
    L = grid.shape[0]
    time_grid = np.zeros((L, L), dtype=int)
    
    # store state fraction at a timestep
    data = {
        'time': [],
        S: [],
        E: [],
        I: [],
        Q: [],
        R: []
    }
    
    # main loop
    for step in range(max_steps):
        # read current timestep
        data['time'].append(step)
        
        # calculate and read current state fraction
        proportions = calculate_state_proportions(grid)
        for state in [S, E, I, Q, R]:
            data[state].append(proportions[state])
        
        # update grid
        grid, time_grid = update_grid(grid, time_grid, params)
    
    return data

# scan any parameter
def parameter_scan(param_name, param_values, base_params, 
                  grid_size=100, population_density=0.75, 
                  init_exposure=50, init_infection=5,
                  max_steps=100, num_runs=1):
   
    import numpy as np
    import time
    
    # validify parameters
    if param_name not in base_params:
        raise ValueError(f"parameter '{param_name}' is not in the base_params")
    
    # create a result dictionary
    results = {}
    
    # store metadata
    metadata = {
        'param_name': param_name,
        'param_values': param_values,
        'base_params': base_params.copy(),
        'grid_size': grid_size,
        'population_density': population_density,
        'init_exposure': init_exposure,
        'init_infection': init_infection,
        'max_steps': max_steps,
        'num_runs': num_runs,
        'run_time': 0,
        'max_infected': {},  # store the max infection rate for each state
        'max_infected_time': {}  # store the timestep that max infection rate appears
    }
    
    print(f"scanning parameter '{param_name}'...")
    start_time_total = time.time()
    
    # simulate every parameter values
    for value in param_values:
        print(f"  simulating {param_name} = {value}")
        # upadate parameter
        params = base_params.copy()
        params[param_name] = value
        
        # init accumulated data
        accumulated_data = {
            'time': list(range(max_steps)),
            S: np.zeros(max_steps),
            E: np.zeros(max_steps),
            I: np.zeros(max_steps),
            Q: np.zeros(max_steps),
            R: np.zeros(max_steps)
        }
        
        # run multiple times and take average
        for run in range(num_runs):
            if num_runs > 1:
                print(f"    run {run+1}/{num_runs}")
                
            # init grid
            grid = initialised_grid(grid_size, population_density, init_exposure, init_infection)
            
            # run simulation
            run_start_time = time.time()
            data = run_simulation(grid, params, max_steps)
            run_end_time = time.time()
            
            # accumulate data
            for state in [S, E, I, Q, R]:
                accumulated_data[state] += np.array(data[state])
        
        # calculate average
        if num_runs > 1:
            for state in [S, E, I, Q, R]:
                accumulated_data[state] /= num_runs
        
        # store result
        results[value] = accumulated_data
        
        # calculate and store max infection rate and time
        max_infected_idx = np.argmax(accumulated_data[I])
        metadata['max_infected'][value] = accumulated_data[I][max_infected_idx]
        metadata['max_infected_time'][value] = accumulated_data['time'][max_infected_idx]
        
        print(f"    max infection ratio: {metadata['max_infected'][value]:.4f} at day {metadata['max_infected_time'][value]}")
    
    # calculate total running time
    end_time_total = time.time()
    metadata['run_time'] = end_time_total - start_time_total
    print(f"finised paramters scanning, time taken: {metadata['run_time']:.2f} seconds")
    
    return results, metadata