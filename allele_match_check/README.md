## Running & Dependencies
To run this project you will need these libraries:

1. openpyxl
2. jupyter notebook (jupyter)
3. matplotlib

Install with pip, for windows users you could use (working for python 3.5+, be sure that py command recognized, if not you will need to change PATH variable):

```
py -m ensurepip
py -m pip install --upgrade pip
py -m pip install openpyxl
py -m pip install jupyter
py -m pip install matplotlib
```

Be sure that project's genetic files configured in main.py, you will need 4 files:

1. The genetic data file, currently supported is excel file format
2. The ambiguities file
3. The haplotypes file
4. The frequency file each line corresponding to the line in the haplotype file that has the same number

Run this project with **python 3.6** by running main.py: 

```py main.py```

Also PyPy recommended (it's reducing running time, install here: https://pypy.org/download.html).

After that you can see the results (solution.csv file by default) using jupyter notebook (there is a notebook file inside notebook folder) and after some configuration.