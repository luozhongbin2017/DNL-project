# Dynamic network loading package for MATLAB
*algorithm by Ke Han | implementation by Gabriel Eve*

### Table of contents
* [Introduction](#introduction)
* [Pre-processor](#pre-processor)
* [Main program](#main-program)
* [GUI](#graphical-display)
* [References](#references)


## Introduction
The package is an implementation of the dynamic network loading (DNL) model formulated as a system of differential algebraic equations proposed by (Han et al., 2014a). This model employs the classical LWR (Lighthill-Whitham-Richards) link dynamics and captures vehicle spillback.

[[back](#table-of-contents "Go to table of contents")]

## Pre-processor
The `create_network_properties.m` script loads a MAT file containing arrays `linkData` and `pathList` , and optional variables: `nodeCoordinates`, `networkName` and `fileName`. The required format of these inputs is described in each section below.

### `linkData`
Array of dimension [*m* x 5], where *m* is number of links in network. Contains information about each link (in rows).

|          | Tail node | Head node | Capacity [veh/s] | Length [m] | Free flow time [s] |
|----------|-----------|-----------|------------------|------------|--------------------|
| link 1   |           |           |                  |            |                    |
| link 2   |           |           |                  |            |                    |
| link 3   |           |           |                  |            |                    |
|          |           |           |                  |            |                    |
| link *m* |           |           |                  |            |                    |

### `pathList`
Array of dimension [*p* x *k*_max] where *p* is the number of paths and *k* is the number of links in each path. *k*_max is the number of links in the longest path, and all paths with *k* < *k*_max must have trailing zeros. Examples shown below (path 1 has 2 links, path 2 has 3 links, and path 3 has *k*_max links):

|          | 1st link | 2nd link | 3rd link | 4th link |   |   | *k*_max link|
|----------|:--------:|:--------:|:--------:|:--------:|:-:|:-:|:-----------:|
| path 1   |     4    |    3     |    0     |    0     |   |   |      0      |
| path 2   |     1    |    2     |    8     |    0     |   |   |      0      |
| path 3   |     8    |    5     |    4     |    7     |   |   |      2      |
|          |          |          |          |          |   |   |             |
| path *p* |     1    |    0     |    0     |    0     |   |   |      0      |


### Optional inputs
`networkName`: name of the network

`fileName`: name to use for file outputs (must meet MATLAB's naming requirements)

`nodeCoordinates`:

|          | X | Y |
|----------|---|---|
| node 1   |   |   |
| node 2   |   |   |
| node 3   |   |   |
|          |   |   |
| node *n* |   |   |



[[back](#table-of-contents "Go to table of contents")]

## Main program

[[back](#table-of-contents "Go to table of contents")]

## Graphical display

[[back](#table-of-contents "Go to table of contents")]

## Versions

[[back](#table-of-contents "Go to table of contents")]

## References

Han, K., Friesz, T.L., Yao, T., 2014a. Vehicle spillback on dynamic traffic networks and what it means for dynamic traffic assignment models. 5th International Symposium on Dynamic Traffic Assignment. Salerno, Italy, 17-19 June 2014.
