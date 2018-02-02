# t4iss
Implementation of some theoretical methods for ISS at NSLS-II.

    # Assuming you have python 3.6+ environment through anaconda 
    $ conda install --channel matsci pymatgen
    $ conda install -c omnia pybtex 
    $ cd data; unzip XANES.zip; cd ..
    $ jupyter notebook &

You need to know your Materials Project API key to run these.

![](img/api.png) 

You should also be familiar with Materials Project ID of a structure.

![](img/mpid.png)



## module-1
This module gets structure from Materials Project based on mpid and generates 
a plot of coordination number around central atom and x-ray absorption spectrum
for each non-equivalent atomic sites.

```python
get_XANES(mpr,mpid='mp-5229',absorbing_atom='Ti',export_figure=True)
```

will generate this plot for Ti-K edge of SrTiO3 (mp-5229):

![](img/mp-5229_Ti.png)

If you want O-K edge, change Ti as O

```python
get_XANES(mpr,mpid='mp-5229',absorbing_atom='O',export_figure=True)
```

You can do alse search in MP database like this:
    
```python
# search in MP
mpid_list = search_MP(mpr,search_pattern='Li-Ti-O',nmax=20)

# plot first 5
for s in mpid_list[0:5]:
    get_XANES(mpr,mpid=s,absorbing_atom='Ti')
```

This will search for string "Li-Ti-O" in MP and retrive ids of structures less than 20 atoms in unit-cell.
First 5 structure will be plotted. Change accordingly....


## Author
* Mehmet Topsakal (mtopsakal@bnl.gov)
