<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# CHANGELOG

## 0.1.2 (2025-04-11)

### Fix

- **ConnectionTable**: show correct values in attributes

### Refactor

- **ConnectionTable**: stop holding HB info after usage

## 0.1.1 (2025-04-11)

### Fix

- resolve numpy RuntimeWarning
- correctly encodes connection information
- correctly open temporary file in MolClusters
- update import path for vdwradii
- resolve type hints errors

## 0.1.0 (2025-04-11)

### Feat

- add attributes to connections

## 0.1.0-dev.0 (2025-04-11)

### Feat

- add radii information to topology
- nucleus analysis
- add clusters nucleus analysis
- add resids properties to Cluster
- add software version and absolute paths to json output
- change dipole to debye
- add json export for molclusters
- add version option to command
- add properties analysis of clusters
- implements initial HB rule algorithm
- initial solute-solvent analysis algorithm
- use YAML file as settings for analysis
- add MolCluster Class and dominance algorithm
- add init files so test can import MolClusters
- finish implementation of ConnTable
- change code to support ConnTable as class
- add generator of sub graphs of ConnTable
- change conn_tab to class ConnTable
- add cluster size evolution analysis
- read large chunks of trajectory at once
- add size property to Cluster
- add trajectory analysis

### Fix

- add encode_connections method to MolClustersData for connection encoding
- avoid error when input does not have solute information
- solve clusters of one molecule
- remove unused variables
- unify use of coord instead of trajectory
- add guards to eliminate errors in evolution calculations
- correct MolClusters to use instance universe
- free command line interface
- import Counter
- consider the reciprocal in input file
- ConnTab raises error if access key do not exist
- change handling of zero distance
- save evolution data to file
- correct position of data in clusters_evo_info
- correct merging of multiple clusters_index
- correctly update conn_tab during trajectory
- correct identification of new clusters
- solve dict size change runtime error
- solve merge and separation problem
- correct identification of clusters

### Refactor

- avoid unnecessary nucleus checks
- move cluster ag and uni to holder class
- break up update_cluster in smaller functions
- modify Cluster to use SubConnTables
- add SubConnTable to ConnTable
- changes ConnTab construction
