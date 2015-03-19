# Dynamic network loading package for MATLAB
*algorithm by Ke Han | implementation by Gabriel Eve*

## Contents
- [Introduction](#introduction)
- [Pre-processor](#pre-processor)
- [Main program](#main-program)
- [GUI](#graphical-display)


## Introduction
An implementation of a dynamic network loading model using the algorithm set out in (missing paper reference) (Han et al., 2012)

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
`nodeCoordinates`:

|          | X | Y |
|----------|---|---|
| node 1   |   |   |
| node 2   |   |   |
| node 3   |   |   |
|          |   |   |
| node *n* |   |   |

`networkName`: name of the network, useful for graphic display of results

`fileName`: name to use for file outputs (must meet MATLAB's naming requirements)

## Main program


## Graphical display

## Versions


