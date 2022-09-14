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
> Installation might take some time since it is compiling ViennaRNA from source.
> 
> Sit back and relax

## Usage

Afterwards you should be able to use the functionality as follows:

```python
from CPExpectedDistance import expected_distance

seq = "AAUGCUGAACGUAGCUACGAUCGAUCGACG"
array = expected_distance(seq)
terminal_distance = array[0, -1]
```

The `expected_distance` function takes a second argument which is a dictionary of 
string keys and values. This is used to change ViennaRNAs model details.
The key names and their value types can be found here:

https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__model__details.html#structvrna__md__s

You can change the folding temperature as follows:

```python
from CPExpectedDistance import expected_distance

seq = "AAUGCUGAACGUAGCUACGAUCGAUCGACG"
array = expected_distance(seq, {"temperature": 20})
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
