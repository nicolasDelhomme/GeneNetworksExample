# Gene Networks Example

## Setup

### Git
```{bash}
git clone https://github.com/nicolasDelhomme/GeneNetworksExample.git
cd GeneNetworksExample
git submodule init
git submodule update --remote
```

### Directories
```{bash}
ln -s ../persistent analysis
```

## Documentation
The `doc` folder contains the introductory lecture and the tutorial.

## Methods
1. The `src/R/dataRetrieval.R` script is used to retrieve the raw data
2. The `src/R/dataPreparation.R` script is used to access the quality of the data, normalise it and export it for the gene network inference.
3. The tutorial in the `doc` folder contains the remaining methods