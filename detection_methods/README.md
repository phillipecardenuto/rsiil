# Detection Methods

This directory contains two CMFD popular libraries that include a feature that modified the detector's output to assign for each detected tampering (object and its copy ) a unique ID, allow us to pinpoint different detection within a image.

## 1. Evaluation of Popular Copy-Move Forgery Detection Approaches [1]
To make a easier installation of this library, we provided a [Dockerfile](vole/Dockerfile) to build this library along with the extra feature



## 2. Automatic detection of internal copy-move forgeries in images [3]

This library is an implementation of the method presented by Cozzolino et al. [2]. 

We also provided a [Dockerfile](patch_match/Dockerfile), aiming a easier installation of this library with the extra feature.

The output of this method is a binary int32 matrix, which the shape is the same as its input, and each suspicious pixels has value different from zero.



## Build all Dockerfiles

```bash 
$ ./build.sh
```



The docker images were also released in Dockerhub:

Dockerhub of the first method [1]:

```
docker pull phillipecardenuto/cmfd:vole
```

Dockerhub of the second method [2,3]:

```
docker pull phillipecardenuto/cmfd:fd
```



## Detection

Check the [Detection Notebook](https://github.com/phillipecardenuto/rsiil/blob/main/notebooks/Detection.ipynb) to see same examples of detection using these dockers.



### BusterNet [4]

In addition to the the previous methods, we also used the [pre-trained](https://github.com/isi-vista/BusterNet) model from [4], during the experiments presented in the paper.



### References

[1] "Evaluation of Popular Copy-Move Forgery Detection Approaches" by V. Christlein, C. Riess, J. Jordan, C. Riess and E. Angelopoulou, IEEE TIFS, vol.
7, no. 6, pp. 1841-1854, Dec. 2012.

[2] Davide Cozzolino, Giovanni Poggi, and Luisa Verdoliva, Efficient dense-field copy-move forgery detection, IEEE Transactions on Information Forensics and Security, 10 (2015),pp. 2284-2297

[3] T. Ehret, “Automatic detection of internal copy-move forgeries in images,” Image Processing On Line, vol. 8, pp. 167–191, Jul. 2018. [Online].
Available: https://doi.org/10.5201/ipol.2018.213

[4] Y. Wu, W. AbdAlmageed, and P. Natarajan, “Busternet: Detecting image copy-move forgery with source/target localization,” in European Conference on Computer Vision (ECCV). Springer, 2018.

