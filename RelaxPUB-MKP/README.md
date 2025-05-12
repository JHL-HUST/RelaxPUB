# RelaxPUB-MKP

The algorithm is implemented on top of DiseMKP (Jiang et al. A Refined Upper Bound and Inprocessing for the Maximum k-plex Problem. IJCAI 2023) with our RelaxPUB.

## Compile the code

```sh
$ ./build
```
It generates an executable "RelaxPUB-MKP".

## Run the code

```sh
$ ./RelaxPUB-MKP {path_to_binary_compressed_graph} -x {k_value}
```

An example of computing the exact maximum 2-plex for the graph brock200-3 is as follows
```sh
$ ./RelaxPUB-MKP ../dataset/Text/brock200-3.col -x 2
```

## Data format
We adopt the text format as DiseMKP did. Several demonstration graphs are given in "../dataset/Text" folder.
