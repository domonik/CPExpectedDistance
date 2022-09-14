# CPExpectedDistance

The CP(ython)ExpectedDistance is a Module to calculate Expected Distances of terminal on the ensemble of RNA Structures.
Under the hood it uses the ViennaRNA C API as well as the numpy API to create an expected Distance matrix.

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

