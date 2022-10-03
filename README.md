# CPExpectedDistance

The CP(ython)ExpectedDistance is a Module to calculate Expected Distances of terminal on the ensemble of RNA Structures.
Under the hood it uses the ViennaRNA C API as well as the numpy API to create an expected Distance matrix.

The idea is from the paper:

```
@article{Clote2012ExpectedDB,
  title={Expected distance between terminal nucleotides of RNA secondary structures},
  author={Peter Clote and Yannick Ponty and Jean-Marc Steyaert},
  journal={Journal of Mathematical Biology},
  year={2012},
  volume={65},
  pages={581-599}
}
```



## Installation

Installation should be as simple as 

```shell
pip install cpexpecteddistance
```

> **Warning**
>   
> Please Note that you need ViennaRNA to be installed. We recommend to install it via Conda
> 
> ```
> conda install -c bioconda viennarna
> ```
> 


> **Warning**
>   
> If you get any compiler errors try to use a conda environment and install g++ there via:
> 
> ```
> conda install -c conda-forge gxx
> ```


## Usage

Afterwards you should be able to use the functionality as follows:

```python
from CPExpectedDistance import expected_distance
import RNA

seq = "AAUGCUGAACGUAGCUACGAUCGAUCGACG"
fc = RNA.fold_compound(seq)
array = expected_distance(fc)
terminal_distance = array[0, -1]
```

The `expected_distance` function is fully compatible with the ViennRNA Python API, thus you can
change the model details for folding as usual. 

For example you can change the folding temperature as follows:

```python
from CPExpectedDistance import expected_distance
import RNA

seq = "AAUGCUGAACGUAGCUACGAUCGAUCGACG"
md = RNA.md(temperature=20)
fc = RNA.fold_compound(seq, md)
array = expected_distance(fc)
terminal_distance = array[0, -1]
```


In case you are using This software for your publications you may want to cite:

```
@article{lorenz2011viennarna,
  title={ViennaRNA Package 2.0},
  author={Lorenz, Ronny and Bernhart, Stephan H and H{\"o}ner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F and Hofacker, Ivo L},
  journal={Algorithms for molecular biology},
  volume={6},
  number={1},
  pages={1--14},
  year={2011},
  publisher={Springer}
}
```
