# Recod.ai Scientific Image Integrity Library

A collection of algorithms to synthetically create scientific image for forensics and integrity analysis.



## Library

The library  implements the most type image tampering functions.

1. Image Duplication
2. Retouching 
3. Cleaning

This [notebook](https://github.com/phillipecardenuto/rsiil/blob/repo-org/notebooks/Tampering-Simple-Scientific-Figures.ipynb) explains how to apply each type of forgery in a scientific image



The library also mimics the behavior of images placed in scientific documents, such as compound figures -- with indicative letters and graphs.

There is two possible types of forgeries for compound figures:

1. Intra-panel (forgeries that are isolated within a single panel from the compound figure)

   [Notebook](https://github.com/phillipecardenuto/rsiil/blob/repo-org/notebooks/Tampering-Compound-Intra-Panel-Scientific-Figures.ipynb) explaining each type of implemented forgery

2. Inter-panel (forgeries that involve more than one panel of the figure):

   [Notebook](https://github.com/phillipecardenuto/rsiil/blob/repo-org/notebooks/Tampering-Compound-Inter-Panel-Scientific-Figures.ipynb) explaining each type of implemented forgery



**Requirements:**

To run the notebooks, make sure to install python3.8 and the modules from included in the [requirements.txt](https://github.com/phillipecardenuto/rsiil/blob/repo-org/requirements.txt) .



# Recod.ai Scientific Image Integrity Dataset

Using the implemented library, we created a synthetic dataset dedicated for forensic purpose of scientific images.

![rsiid](https://github.com/phillipecardenuto/rsiil/blob/repo-org/.figs/rsiid.jpg)

*Dataset*:

[Train set](https://drive.google.com/file/d/1ueRs8ySiFCFEcsIiB3X6KiQcd2aVBMbj/view?usp=sharing)

[Test set](https://drive.google.com/file/d/1E3HHXSbjIMnfiFK9XxSvSRa6RcUw2-SD/view?usp=sharing)



*Source figures* and *compound figure templates* used to create the tampering dataset:

[Source Figures]( https://drive.google.com/file/d/1l5RdipKFm_bnDxSwRYwJVG2NUYij7T3w/view?usp=sharinG)

[Templates](https://drive.google.com/file/d/11DXW_yFEZsje0aJyZDGaT4cJwwz2Dmzt/view?usp=sharinG)



## Dataset Organization

Both train and test set have simple and compound figures, organized with the following schematic:

**Simple Images**

![](https://github.com/phillipecardenuto/rsiil/blob/repo-org/.figs/simple-data.jpg)

**Compound figure**

![](https://github.com/phillipecardenuto/rsiil/blob/repo-org/.figs/compound-data.jpg)





### Citation

The dataset is distributed under [Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.en)

If you use any of this content for research purpose, please cite:

```tex
 @article{cardenuto_2022, 
 title={Benchmarking scientific image forgery detectors},
 volume={28}, DOI={10.1007/s11948-022-00391-4},
 number={4},
 journal={Science and Engineering Ethics},
 author={Cardenuto, João P. and Rocha, Anderson}, year={2022}
 } 
```

