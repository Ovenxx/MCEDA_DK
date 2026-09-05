# Multi-Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi-Objective Green Vehicle Routing

> This repository is the official code release for the paper: **W. Wang, Y. Shen, Y. Li, and Y. Fan, "Multi Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi Objective Green Vehicle Routing," *IEEE Transactions on Evolutionary Computation*, 2026 (Early Access). [DOI: 10.1109/TEVC.2026.3704658](https://doi.org/10.1109/TEVC.2026.3704658)**

## Overview

This paper formulates a multi-objective time-space-dependent green VRP with pickup-delivery and time window (MO-TSD-GVRPPDTW), which simultaneously incorporates pickup demands, delivery demands, and time-space dependences. To solve this problem, a multi-consensus-based estimation of distribution algorithm embedded with domain knowledge (MCEDA-DK) is proposed.

It comprises three key components: 1) a domain knowledge embedding (DKE) method; 2) a multi-consensus evolution (MCE) strategy; and 3) a multi-scale neighborhood search (MNS) strategy.

The DKE method constructs a time-space compatibility matrix to prune the search space, improving both initialization quality and evolution efficiency.

The MCE strategy iteratively updates the population's consensus on multiple conflicting objectives to indicate diversified non-dominated combinations, guiding subsequent searches.

The MNS strategy stimulates the population to search both large- and small-scale neighborhoods, enabling persistent exploration of better local regions and identification of local non-dominated solutions, which effectively balances diversified and intensified searches.

## Quick Start: Reproducing the Results

### 1. Clone the repository
```bash
git clone https://github.com/Ovenxx/Multi-objective-GVRP.git
cd Multi-objective-GVRP
```

### 2. Reproduce benchmark comparison results
1) Open MATLAB and set the current folder to the Benchmark experiments/ directory.
```matlab
>> cd('Benchmark experiments')
```

2) Run the experiment script to generate the non-dominated solution sets:
```matlab
>> Experiment_script_100
```

3) Compute the hypervolume (HV) values and obtain the rHV variable (Its calculation employs Monte Carlo sampling):
```matlab
>> HV_comparison
```

4) Print the comparison result:
```matlab
>> fprintf('AVERAGED HV VALUES OBTAINED BY MCEDA-DK, MOLS, AND MOMA ON THE WC INSTANCE SET.\n');
>> disp(rHV);
```

## Citation

If you use this code or find this repository helpful in your research, please cite the following paper:

> W. Wang, Y. Shen, Y. Li, and Y. Fan, "Multi Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi Objective Green Vehicle Routing," *IEEE Transactions on Evolutionary Computation*, 2026 (Early Access), doi: [10.1109/TEVC.2026.3704658](https://doi.org/10.1109/TEVC.2026.3704658).

**BibTeX:**

```bibtex
@article{wang2026multi,
  author  = {Wang, Wei and Shen, Yindong and Li, Yang and Fan, Yibao},
  title   = {Multi Consensus Based Estimation of Distribution Algorithm with Domain Knowledge for Multi Objective Green Vehicle Routing},
  journal = {IEEE Transactions on Evolutionary Computation},
  year    = {2026},
  doi     = {10.1109/TEVC.2026.3704658},
  note    = {Early access}
}
```

A [CITATION.cff](CITATION.cff) file is also provided, so GitHub's "Cite this repository" feature works out of the box.
