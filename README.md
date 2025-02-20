# Recod.ai Scientific Image Integrity Library

A collection of algorithms to synthetically create scientific images for forensics and integrity analysis.

An in-depth explanation of each algorithm and dataset is described in our research work:
[Benchmarking Scientific Image Forgery Detectors](https://link.springer.com/article/10.1007/s11948-022-00391-4)



## Library

The library implements the most type of image tampering functions.

1. Image Duplication
2. Retouching 
3. Cleaning

This [notebook](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Simple-Scientific-Figures.ipynb) explains how to apply each type of forgery in a scientific image.


The library also mimics the behavior of images placed in scientific documents, such as compound figures -- with indicative letters and graphs.

There are two possible types of forgeries for compound figures:

1. Intra-panel (forgeries that are isolated within a single panel from the compound figure)

   [Notebook](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Compound-Intra-Panel-Scientific-Figures.ipynb) explaining each type of implemented forgery

2. Inter-panel (forgeries that involve more than one figure panel):

   [Notebook](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Compound-Inter-Panel-Scientific-Figures.ipynb) explaining each type of implemented forgery



**Requirements:**

To run the notebooks, make sure to install python3.8 and the modules included in the [requirements.txt](https://github.com/phillipecardenuto/rsiil/blob/main/requirements.txt).



# Recod.ai Scientific Image Integrity Dataset

Using the implemented library, we created a synthetic dataset dedicated to forensics purposes and scientific integrity.

![rsiid](https://github.com/phillipecardenuto/rsiil/blob/main/.figs/rsiid.jpg)

*[RSIID](http://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/) Dataset*:

[Train set](http://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/trainset.zip)

[Test set](http://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/testset.zip)

If you are facing issues while downloading or accessing these files, try:
```bash
curl -k https://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/trainset.zip --output trainset.zip
curl -k https://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/testset.zip --output testset.zip
```



*Source figures* and *compound figure templates* used to create the tampering dataset:

[Source Figures](http://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/artificial_forgery_src_data.zip)

[Templates](http://intranet.recod.ic.unicamp.br/~jcardenuto/pub/rsiid/template.zip)



## Dataset Organization

Both train and test sets have simple and compound figures, organized with the following schematic:

**Simple Images**

![](https://github.com/phillipecardenuto/rsiil/blob/main/.figs/simple-data.jpg)

**Compound figure**

![](https://github.com/phillipecardenuto/rsiil/blob/main/.figs/compound-data.jpg)





### Citation

The dataset is distributed under the [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.en) license.

If you use any content from this repository, please cite:

```tex
 @article{cardenuto_2022, 
 title={Benchmarking scientific image forgery detectors},
 volume={28}, DOI={10.1007/s11948-022-00391-4},
 number={4},
 journal={Science and Engineering Ethics},
 author={Cardenuto, João P. and Rocha, Anderson}, year={2022}
 } 
```

