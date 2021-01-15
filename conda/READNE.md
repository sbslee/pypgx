```
cd pypgx
cd conda
conda-build pypgx .
```

```
conda convert --platform linux-64 /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2
```

```
anaconda upload /Users/sbslee/opt/anaconda3/conda-bld/osx-64/pypgx-0.1.37-py38_0.tar.bz2
anaconda upload linux-64/pypgx-0.1.37-py38_0.tar.bz2

- requests
- pandas
- bs4
- lxml
- pysam
- jupyter_core
```
