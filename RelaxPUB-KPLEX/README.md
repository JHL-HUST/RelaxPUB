# RelaxPUB-KPLEX

The algorithm is implemented on top of KPLEX (Wang et al. A Fast Maximum k-Plex Algorithm Parameterized by the Degeneracy Gap. IJCAI 2023) with our RelaxPUB.

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "RelaxPUB-KPLEX".

## Run the code

```sh
$ ./RelaxPUB-KPLEX {path_to_binary_compressed_graph} {k_value}
```

An example of computing the exact maximum 2-plex for the graph brock200-3 is as follows
```sh
$ ./RelaxPUB-KPLEX ../dataset/Binary/brock200-3.bin 2
```

The solution is in output file "kplexes.txt".

## Data format
We adopt the time-efficient binary format rather than the text format.  Several demonstration graphs are given in "../dataset/Binary" folder.

Transforming [network-repo graphs](http://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz) from text format to binary format is as follows:
```sh
$ g++ ./toBin.cpp -o toBin
$ mv ../dataset/Text/brock200-3.col ./
$ ./toBin brock200-3.col
$ ls brock200-3.bin
```