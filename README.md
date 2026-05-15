This repository contains the official implementation of the paper:
**"Multi-Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi-Objective Green Vehicle Routing"**  

## Overview
This paper formulates a multi-objective time-space-dependent green VRP with pickup-delivery and time window (MO-TSD-GVRPPDTW), which simultaneously incorporates pickup demands, delivery demands, and time-space dependences. To solve this problem, a multi-consensus-based estimation of distribution algorithm embedded with domain knowledge (MCEDA-DK) is proposed. It comprises three key components: 1) a domain knowledge embedding (DKE) method; 2) a multi-consensus evolution (MCE) strategy; and 3) a multi-scale neighborhood search (MNS) strategy. The DKE method constructs a time-space compatibility matrix to prune the search space, improving both initialization quality and evolution efficiency. The MCE strategy iteratively updates the population’s consensus on multiple conflicting objectives to indicate diversified non-dominated combinations, guiding subsequent searches. The MNS strategy stimulates the population to search both large- and small-scale neighborhoods, enabling persistent exploration of better local regions and identification of local non-dominated solutions, which effectively balances diversified and intensified searches.

## Repository Structure
.
├── Benchmark experiments/ # Experiments on Solomon benchmark instances
│ ├── MCEDA_DK.m # Core MCEDA-DK algorithm implementation
│ ├── Experiment_script_100.m # Main script for 100-customer instances
│ ├── HV_comparison.m # Hypervolume computation and comparison
│ ├── hypervolume.m # Hypervolume calculation utility
│ ├── cdp101_data.mat # Problem instance data (100 customers)
│ ├── cdp201_data.mat
│ ├── benchmark_Solomon/ # Solomon benchmark datasets
│ ├── MCEDA_DK_data/ # Results of MCEDA-DK
│ ├── MD_MOLS_data/ # Results of comparison algorithm MD-MOLS
│ └── MD_MOMA_data/ # Results of comparison algorithm MD-MOMA
└── JD Logistics experiments/ # Real-world JD Logistics case studies
├── Experiment_script_200.m # Script for 200-customer instances
├── Experiment_script_400.m # Script for 400-customer instances
├── JD_2001_data.mat # Logistics data (200 customers)
├── JD_2002_data.mat
├── JD_2003_data.mat
├── JD_2004_data.mat
├── JD_4001_data.mat # Logistics data (400 customers)
├── JD_4002_data.mat
└── Liu_Tang_Yao/ # JD Logistics datasets

## Quick Start: Reproducing the Results

### 1. Clone the repository
```bash
git clone https://github.com/Ovenxx/MCEDA_DK.git
cd MCEDA_DK

### 2. Reproduce benchmark comparison results
1) Open MATLAB and set the Current Folder to the Benchmark experiments/ directory inside the cloned repository.

2) Run the experiment script to generate the non-dominated solution sets for MCEDA-DK:
```matlab
>> Experiment_script_100

3) Compute the hypervolume (HV) values and obtain the rHV variable:
```matlab
>> HV_comparison

4) Print the comparison result:
```matlab
>> fprintf('AVERAGED HV VALUES OBTAINED BY MCEDA-DK, MOLS, AND MOMA ON THE WC INSTANCE SET.\n');
>> disp(rHV);

