# RelaxGCB-MKP

The algorithm is implemented on top of DiseMKP (Jiang et al. A Refined Upper Bound and Inprocessing for the Maximum k-plex Problem. IJCAI 2023) with our RelaxGCB.

## Compile the code

```sh
$ ./build
```
It generates an executable "RelaxGCB-MKP".

## Run the code

```sh
$ ./RelaxGCB-MKP {path_to_binary_compressed_graph} -x {k_value}
```

An example of computing the exact maximum 2-plex for the graph brock200-3 is as follows
```sh
$ ./RelaxGCB-MKP ../dataset/Text/brock200-3.col -x 2
```

## Data format
We adopt the text format as DiseMKP did. Several demonstration graphs are given in "../dataset/Text" folder.
