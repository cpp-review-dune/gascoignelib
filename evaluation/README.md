# Evaluation
In this folder you can find a evaluation and testing skript. 
It will build all the examples in the 'examples' dir and run them with increasing refinement levels of the multigrid.

To run the test run:
```
python3 -m unittest evaluation/evaluation.py
```
in the root dir of this repository.

For a plotting example run:
```
python3 evaluation/evaluation.py
```
This will run all the examples and plot a timer example, see the end of the evaluation.py.