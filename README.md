# Using an SEIQR model to show how the spread of a viral infection affects a population

## Background
  This Project is originally from the Module NSCI0036: Programming for Scientists 24/25 at department of Natural Sciences, UCL. Now the Project outcomes including Python codes and project report are released to you to practise if you are interested in Celluar Automation Model writing in Python.

## Introduction
  The early days of COVID-19 demonstrated how rapidly viruses can spread and their impact on population density. The first groupof patients were identified on the 21st December 2019 in the city of Wuhan, which had a population of 11 million people. By 25thJanuary 2020, an estimated 75,815 individuals had been infected, exhibiting just how quickly an epidemic can spread, escalate,and cause a global health crisis (Qu et al., 2020). In this investigation, we’ll explore the spread of COVID-19 through apopulation, to better understand its epidemiology and spread using the complex computational model SEIQR in Python. Thisstimulates, visualises and forecasts the spread of an epidemic and its effect on population density. This representation of spread isbased on the more simplistic, probabilistic SIR model within a cellular automaton (CA) and demonstrates how populations canchange states corresponding to the chain of infection where the spread is proliferated longitudinally. SIR only employs the use ofthree states: susceptible, infected and recovered, whereas the SEIQR model includes additional states to create a variant modelwhich promises a stronger precision in simulation. Therefore, the SEIQR is the model better suited to investigate the effect ofspread of a viral infection on population density and effect of lockdown rules in different stages. This report will focus on how theinteraction range affects viral spread and how the timing of regional lockdown measures impacts the pattern of transmission,where both the intrinsic nature of infection spread, and governmental epidemic management are considered.

  Please see the references in the report.

  ## Method
  ### SEIQR Model Description

The SEIQR model is initialised by defining the grid parameters, including grid size, population density, and the initial number of exposed and infected individuals. The model is implemented as a cell lattice of size $L \times L$, which in this context represents the population size. Each cell corresponds to an individual and is assigned one of five states, representing distinct stages in the infection chain. The state of each cell evolves over time based on interactions with its neighbouring cells (Table 1).

---

### State Definitions and Transition Probabilities

| State | Description | Transition Probability |
|------|------------|------------------------|
| Susceptible ($S$) | Uninfected, can become exposed | $P_e = 0.5$ |
| Exposed ($E$) | In contact with infected, not yet infectious | $P_i = 0.5$ after $T_i = 7$ days |
| Infected ($I$) | Actively infectious | $P_q = 0.1$ for quarantine |
| Quarantined ($Q$) | Isolated, still infected | $P_r = 0.12$ for recovery |
| Removed ($R$) | Recovered with immunity, or passed away | After $T_r = 21$ days |

**Table 1:** SEIQR model states, their meanings, and associated probability values.

---

Each cell interacts with a fixed number of neighbouring cells defined by the **Moore neighbourhood**. As timesteps progress, each susceptible cell with at least one exposed or infected neighbour may become exposed according to the probability $P_e$. If the transition condition is not satisfied, the individual remains susceptible and uninfected.

---

### Model Parameters

| Parameter | Value | Description |
|----------|-------|-------------|
| $L \times L$ | $100 \times 100$ | Grid size, total population $= 10{,}000$ |
| $D$ | $1$ | Moore neighbourhood radius (subject to change) |
| $P_e$ | $0.5$ | Probability of exposure |
| $P_i$ | $0.5$ | Probability of infection |
| $P_q$ | $0.1$ | Probability of quarantine |
| $P_r$ | $0.12$ | Probability of recovery |
| $T_i$ | $7$ days | Incubation period |
| $T_q$ | $7$ days | Time spent in quarantine |
| $T_r$ | $21$ days | Recovery period |

**Table 2:** Model parameters, assigned values, and descriptions.

---

To implement the cellular automaton (CA), all state transition probabilities are defined as numerical values (Table 2). Cell interactions are governed by the Moore neighbourhood with interaction radius $D$. In the code, $P_e$ represents the probability of exposure, $P_i$ the probability of infection, $P_q$ the probability of quarantine, and $P_r$ the probability of recovery.

The system evolves iteratively using loop-based updates, allowing cell states to change according to predefined transition rules (Table 3). The simulation proceeds over discrete timesteps, incorporating an incubation period $T_i$, a quarantine period $T_q$, and a recovery period $T_r$.

---

### Allowed State Transitions

| Transition | Description |
|-----------|------------|
| $S \rightarrow E$ | If an individual has an adjacent exposed or infected neighbour, with probability $P_e$ |
| $E \rightarrow I$ | After $T_i$, an exposed individual becomes infected with probability $P_i$ |
| $I \rightarrow Q$ or $R$ | After $T_i + T_q$, an infected individual may be quarantined with probability $P_q$ |
| $I \rightarrow R$ | After $T_i + T_q$, an infected individual may recover with probability $P_r$ |
| $Q \rightarrow R$ | After $T_r$ timesteps, a quarantined individual recovers with probability $P_r$ |

**Table 3:** Allowed state transitions and corresponding rules.

---

### Simulation Procedure

The simulation proceeds through discrete timesteps, each consisting of the following three steps:

#### 1. Compute state proportions
At the beginning of each timestep, the proportions of individuals in each state ($E$, $I$, $Q$, $R$) are calculated. The number of cells in each state is counted and normalised by the total population to obtain fractional values, which illustrate the progression of the disease.

#### 2. Apply state transition rules
Each cell is evaluated to determine whether it can change state according to the transition rules listed in Table 3. Disease spread is governed by interactions within the Moore neighbourhood and the interaction range $D$ (Davies, 1995). Each susceptible individual assesses its surroundings by counting neighbouring exposed and infected cells.

#### 3. Record data for analysis
A separate grid is used to track the duration each cell remains in a given state. This ensures that the parameters $T_i$, $T_q$, and $T_r$ are enforced correctly, resetting only when a state transition occurs. At the end of each timestep, the updated state proportions are recorded, enabling the generation of plots such as temporal line graphs and spatial simulation visualisations presented in this report.
