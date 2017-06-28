# IDBA to GFA

This is a simple tool for converting [IDBA](https://github.com/loneknightpy/idba) assemblies into [GFA graphs](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md). This allows the graphs to be loaded into [Bandage](https://github.com/rrwick/Bandage) or anything else which accepts GFA files.

It was made my myself and [Matt Olm](https://github.com/MrOlm).


### Requirements

It's just a single Python 3 script, so no installation is required. It uses the [`print_graph`](https://github.com/loneknightpy/idba/blob/master/src/tools/print_graph.cpp) tool (which is built along with IDBA), so that must be installed on the same computer.


### Usage

```
usage: idba_to_gfa.py [-h] [--print_graph PRINT_GRAPH] idba_assembly kmer

IDBA to GFA: a tool for converting IDBA assemblies to GFA graphs

positional arguments:
  idba_assembly         Assembly from IDBA, e.g. scaffold.fa
  kmer                  assembly k-mer size (e.g. 100)

optional arguments:
  -h, --help            show this help message and exit
  --print_graph PRINT_GRAPH
                        Location of IDBA's print_graph tool (required if
                        print_graph is not in PATH)

output: GFA to stdout
```


### Example commands

Build a 100-mer graph from IDBA scaffolds:<br>
`idba_to_gfa.py scaffold.fa 100 > scaffold_100.gfa`

Build an 80-mer graph from IDBA contigs:<br>
`idba_to_gfa.py contig.fa 80 > contig_80.gfa`

Give the location of `print_graph` (necessary if it's not in PATH):<br>
`idba_to_gfa.py --print_graph /path/to/idba-1.1.3/bin/print_graph scaffold.fa 100 > scaffold_100.gfa`


### License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
