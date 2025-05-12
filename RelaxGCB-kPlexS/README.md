# RelaxGCB-kPlexS

The algorithm is implemented on top of kPlexS proposed in the VLDB 2023 paper shown below with our RelaxGCB.

<pre>
Lijun Chang, Mouyi Xu and Darren Strash.
<a href="https://lijunchang.github.io/pdf/2022-Maximum-kPlex.pdf">Efficient Maximum k-Plex Computation over Large Sparse Graphs.</a>
Proc. VLDB Endow. 16(2), (2022)
</pre>

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "RelaxGCB-kPlexS".

## Run the code

```sh
$ ./RelaxGCB-kPlexS -g {path_to_graph} -a exact -k {k_value} -o
```

An example of computing the exact maximum 2-plex for the graph brock200-3 is as follows
```sh
$ ./RelaxGCB-kPlexS -g ../dataset/Edges/brock200-3 -a exact -k 2 -o
```

## Data format
Two data formats are supported. The default data format is "edges.txt", which contains a list of undirected edges represented as vertex pairs. The first line contains two numbers n and m, representing the number of vertices and the number of undirected edges, respectively. Note that, the vertex ids must be between 0 and n-1.
