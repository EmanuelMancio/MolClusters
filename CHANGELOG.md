<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# CHANGELOG

## Unreleased

### Fix

- **release**: keep uv.lock in sync and use GitHub's exact skip-ci token

## v0.6.0 (2026-09-22)

### Feat

- **config**: warn when lammps_resnames ranges overlap for the same name
- **config**: add lammps_resnames mapping for LAMMPS topologies
- also write per-cluster-identity coordinate output

### Fix

- **ci**: pin setup-uv to an exact tag, v10 doesn't exist
- **pre-commit**: scope reuse-lint-file hook to the pre-commit stage
- temp-file hygiene, logging consistency, and TypeVar bound
- **config**: store lammps_resnames as ranges instead of per-id entries
- correct dominance-adjacent correctness bugs from code review
- correct config drift and untrack test fixtures gitignore rule
- restrict reuse-lint-file hook to pre-commit stage

### Refactor

- **config**: simplify rule and solute-keyword parsing

## v0.5.0 (2026-09-22)

### BREAKING CHANGE

- dropped python 3.11 support.

### Feat

- **config**: implement MolClsConfig class to handle configuration
- **log**: allow string to be passed to logger initializer
- **log**: change log decorator so no parentheses are needed
- **log**: add logging tools and configuration
- drop python 3.11 and add new dependencies
- **main**: add trajectory loading options for memory management
- **MolClusters**: add support for ignore_composition in cluster analysis
- **SubConnTable**: add resnames property
- add attributes to connections
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

- correct package name in version
- update config keys from 'solute' to 'nucleus' and adjust JSON indentation
- use upper element name for vdwradii
- **main**: correct element extraction from atom name
- **ConnectionTable**: show correct values in attributes
- resolve numpy RuntimeWarning
- correctly encodes connection information
- correctly open temporary file in MolClusters
- update import path for vdwradii
- resolve type hints errors
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

- update type hints to use built-in types (e.g., list, dict) across multiple files
- **ConnectionTable**: stop holding HB info after usage
- avoid unnecessary nucleus checks
- move cluster ag and uni to holder class
- break up update_cluster in smaller functions
- modify Cluster to use SubConnTables
- add SubConnTable to ConnTable
- changes ConnTab construction

## v0.4.0 (2026-05-09)

### BREAKING CHANGE

- dropped python 3.11 support.

### Feat

- **config**: implement MolClsConfig class to handle configuration
- **log**: allow string to be passed to logger initializer
- **log**: change log decorator so no parentheses are needed
- **log**: add logging tools and configuration
- drop python 3.11 and add new dependencies

### Fix

- correct package name in version

### Refactor

- update type hints to use built-in types (e.g., list, dict) across multiple files

## v0.3.1 (2026-03-03)

### Fix

- update config keys from 'solute' to 'nucleus' and adjust JSON indentation
- use upper element name for vdwradii

## v0.3.0 (2025-04-11)

### Feat

- **main**: add trajectory loading options for memory management

### Fix

- **main**: correct element extraction from atom name

## v0.2.0 (2025-04-11)

### Feat

- **MolClusters**: add support for ignore_composition in cluster analysis
- **SubConnTable**: add resnames property

## v0.1.2 (2025-04-11)

### Fix

- **ConnectionTable**: show correct values in attributes

### Refactor

- **ConnectionTable**: stop holding HB info after usage

## v0.1.1 (2025-04-11)

### Fix

- resolve numpy RuntimeWarning
- correctly encodes connection information
- correctly open temporary file in MolClusters
- update import path for vdwradii
- resolve type hints errors

## v0.1.0 (2025-04-11)

### Feat

- add attributes to connections

## v0.1.0-dev.0 (2025-04-11)

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
