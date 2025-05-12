# RelaxGCB-Maplex

The algorithm is implemented on top of Maplex (Zhou et al. Improving Maximum k-plex Solver via Second-Order Reduction and Graph Color Bounding. AAAI 2021) with our RelaxGCB.

## Compile the code

```sh
$ ./build
```
It generates an executable "RelaxGCB-Maplex".

## Run the code

```sh
$ ./RelaxGCB-Maplex {path_to_binary_compressed_graph} {k_value} {cut-off time}
```

An example of computing the exact maximum 2-plex for the graph brock200-3 is as follows
```sh
$ ./RelaxGCB-Maplex ../dataset/Binary/brock200-3.bin 2 1800
```

## Data format
We adopt the time-efficient binary format rather than the text format.  Several demonstration graphs are given in "../dataset/Binary" folder.
