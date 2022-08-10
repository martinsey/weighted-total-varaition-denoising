# Weighted Total Varaition Denoising
A matlab package based on the [paper] (http://www2.mathematik.hu-berlin.de/~rautenca/ContentHP/Prints/Bilevel2_JMIV.pdf) using a staggered grid for the purpose of denoising.

## Usage
Source code is in the `src` folder. For denoising execute
```
cd src
bilevelTV_denoising(f, sigma)
```
where
* f the corrupted image in matlab format. 
* sigma noise level
for example the image can be loaded via can be loaded via
```
load data/input/cam.mat
f = cam_noisy_01
sigma=0.1
```
or for an image from `.png, .jpeg, ...` format
```
f=readimage($PATH_TO_IMAGE)
f=im2double(f);
f=rgb2gray(f);
```
Results will be written in `data/output` as `result.mat` file, which contains the reconstruction and the alpha.

Note that execution for one image can take a while. On a macbook pro 2016 it took 1 hour.

## Structure
The source code will be in the src file, while the actual execution is in the `bilevelTV_denoising.m` file. All data generated will go into
data/output. Some sample input images are in `data/input`. Some pregenerated examples of input and output are in the `data/samples` folder.
