# t4iss
Implementation of some theoretical methods for ISS at NSLS-II.

    $ conda install --channel matsci pymatgen
    $ conda install -c omnia pybtex 
    $ cd data; unzip data.zip; cd ..
    $ jupyter notebook --no-browser &

## module-1
This module gets structure from Materials Project based on mpid and generates 
a plot of coordination number around central atom and x-ray absorption spectrum
for each non-equivalent atomic sites.

![](img/mp.png)


## Author
* Mehmet Topsakal (mtopsakal@bnl.gov)
