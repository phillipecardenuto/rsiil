# Recod.ai Scientific Image Integrity Library (RSIIL)

A collection of algorithms to synthetically create scientific images for forensics and integrity analysis.

## Research Publication

For an in-depth explanation of the algorithms and dataset, please refer to our research paper:
* **Benchmarking Scientific Image Forgery Detectors:** [https://link.springer.com/article/10.1007/s11948-022-00391-4](https://link.springer.com/article/10.1007/s11948-022-00391-4)

## Scientific Image Tampering Library

The library implements various image tampering functions commonly found in scientific image manipulation:

1.  **Image Duplication**
2.  **Retouching**
3.  **Cleaning**

* **Notebook:** [Tampering Simple Scientific Figures](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Simple-Scientific-Figures.ipynb) - Demonstrates how to apply each type of forgery.

The library also supports the creation of compound figures, mimicking the structure of images in scientific documents. Compound figures can have two types of forgeries:

1.  **Intra-panel Forgeries:** Forgeries within a single panel.
    * **Notebook:** [Tampering Compound Intra-Panel Scientific Figures](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Compound-Intra-Panel-Scientific-Figures.ipynb)
2.  **Inter-panel Forgeries:** Forgeries involving multiple panels.
    * **Notebook:** [Tampering Compound Inter-Panel Scientific Figures](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Tampering-Compound-Inter-Panel-Scientific-Figures.ipynb)

**Requirements**

* Python 3.8
* Install the required modules from [requirements.txt](https://github.com/phillipecardenuto/rsiil/blob/main/requirements.txt).

# Recod.ai Scientific Image Integrity Dataset (RSIID)

This dataset, created using the library, is designed for forensics and scientific integrity research.

<div style="text-align: center;">
  <img src="https://github.com/phillipecardenuto/rsiil/blob/main/.figs/rsiid.jpg" width="400">
</div>

* **Dataset Link (Zenodo):** [https://zenodo.org/records/15095089](https://zenodo.org/records/15095089) - Contains all images, metadata, and related files.


## Dataset Organization

The dataset includes both simple and compound figures in the training and testing sets.

**Simple Images**

<img src="https://github.com/phillipecardenuto/rsiil/blob/main/.figs/simple-data.jpg" width="400">


**Compound figure**

<img src="https://github.com/phillipecardenuto/rsiil/blob/main/.figs/compound-data.jpg" width="400">




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

