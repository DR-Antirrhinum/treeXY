# treeXY
## Introduction
treeXY is a software tool for exploring patterns of genetic diversity between groups of taxa mapped to a common reference genome. In its simplest use case, treeXY computes a series of population genetic measures (d<sub>XY</sub>, π<sub>w</sub>, F<sub>ST</sub>, and D) between all taxa in a pairwise manner. These stats are reported in window averages of a user-defined size and overlap.

The novelty of this software is in the computation of UPGMA trees based on the d<sub>XY</sub> distribution at biallelic sites within the genome. These trees facilitate the summary of genetic distance relationships between multiple taxa simultaneously, allowing for the identification of genomic regions showing elevated divergence between groups of taxa, which may underlie genetic barriers. An advantage of d<sub>XY</sub> is that, unlike F<sub>ST</sub>, it is insensitive to changes to π<sub>w</sub> following a population split. This means that it will not be confounded by population effects such as historical selective sweeps, or population bottlenecks.

To summarise d<sub>XY</sub> trees, the following statistics are reported:

| Example                           | Statistic                  | Description                                                                               |
| --------------------------------- |:--------------------------:|:-----------------------------------------------------------------------------------------:|
| ![srb_tree_example](SRB_tree.png) | Shortest Root Branch (SRB) | The distance between the top of the tree, and the higher of the two outermost clades      |
| IMAGE - Root Division example     | Root Division              | A numerical representation of the two outermost clades, when the tree is split at the top |

## Preparing data for analysis
The basic pipeline for preparing data for analysis is as follows:

```mermaid
graph LR;
    FASTQ_1-->BAM_1-->Pileup;
    FASTQ_2-->|BWA| BAM_2-->|Samtools mpileup| Pileup-->|mpileup2sync| SYNC;
    FASTQ_3-->BAM_3-->Pileup;
    
```

The input for treeXY is a Popoolation2 SYNC file. This is a compact representation of the common Pileup format, which contains allele depth information across genomic regions.

## Installation
treeXY requires Python3, along with numpy, pandas, and scipy. It has been tested on the following versions:
- `python-3.9`
- `numpy-1.23.4`
- `pandas-1.5.1`
- `scipy-1.9.3`

A Dockerfile is included for ease of installation. Alternatively, use a Python dependency manager such as Poetry. Use other versions with caution.

## Arguments
`-f` `--file` Input SYNC file to be processed.

`-m` `--min_depth` Minimum depth, across all populations, for a site to be included. The default value is 15.

`-M` `--max_depth` Maximum depth, across all populations, for a site to be included. The default value is 200.

`-a` `--min_allele_depth` Minimum depth for an allele call. The default value is 2.

`-A` `--min_allele_pops` Minimum number of populations required for an allele call. The default value is 2.

`-w` `--window_size` Size of sliding window. The default value is 10000.

`-o` `--window_overlap` Overlap of consecutive sliding windows. The default value is 9000.

`--ignore_multiallelic` Ignore sites with >2 alleles. Normally, treeXY will process these sites by removing the least common allele.

`--write_sync` Write sites passing all treeXY depth checks to a new SYNC file.

`--dxy_trees` If specified, treeXY will run hierarchical clustering to generate dXY trees at all valid biallelic sites. dXY tree statistics, incuding Tree Height, Shortest Root Branch (SRB), and Root Division representation, will be written to a separate CSV file. See Introduction for more details.

`--d_trees` As above, except using Nei's D.

## Usage example
treeXY has been used to study barriers to gene flow between two subspecies of the snapdragon *Antirrhinum majus*, *A. m. pseudomajus* and *A. m. striatum* (Richardson *et al.*, manuscript in preparation).

```geojson 
{ "type": "FeatureCollection",
    "features": [
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [1.772739, 42.763361]},
        "properties": {"population_name": "UNA"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [1.568953, 42.869186]},
        "properties": {"population_name": "BED"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.260464, 42.968486]},
        "properties": {"population_name": "LU"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.2231055, 42.798143]},
        "properties": {"population_name": "AXA"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.039864, 42.725164]},
        "properties": {"population_name": "MIJ"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.122297, 42.507878]},
        "properties": {"population_name": "MON"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.8552415, 42.467675]},
        "properties": {"population_name": "PER"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.58705, 42.643378]},
        "properties": {"population_name": "BOU"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.367453, 42.587006]},
        "properties": {"population_name": "VIL"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.4876195, 42.3895975]},
        "properties": {"population_name": "ARS"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.721694, 42.644139]},
        "properties": {"population_name": "THU"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [3.124183, 42.489458]},
        "properties": {"population_name": "BAN"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.6084845, 42.4479485]},
        "properties": {"population_name": "ARL"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [1.533579, 43.311569]},
        "properties": {"population_name": "CIN"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.052929, 42.326943]},
        "properties": {"population_name": "YP1"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [1.926958, 42.359921]},
        "properties": {"population_name": "YP4"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.091375, 42.322234]},
        "properties": {"population_name": "MP4"}
      },
      { "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [2.170284, 42.331038]},
        "properties": {"population_name": "MP11"}
      }
   ]
}
```

These 18 populations were sampled, sequenced, and mapped to the *Antirrhinum majus* reference genome.

Companions scripts for this analysis can be found in the phylogenetic_forests repository