---
layout: page
title: "Install"
published: true
categories: []
tags: []
---

# Compiling From Source

## Dependencies

Required:

* [CMake](https://cmake.org)
* [Arnold SDK](https://www.solidangle.com/arnold/download/)

Optional:

* [Visual Studio Community](https://www.visualstudio.com/products/visual-studio-community-vs) 2013 or later
* gcc-4.9 or later
* [OpenImageIO](https://github.com/OpenImageIO/oiio) (for image viewing and comparing)

## Windows

1. clone source code from github
<samp>git clone https://github.com/shihchinw/rlShaders.git source-root</samp>
2. open `msvc` command prompt (e.g. _VS2013 x64 Native Tools Command Prompt_) and run `cmake` command for configuration:
<samp>cd source-root
mkdir build
cd build
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DARNOLD_ROOT=path-to-arnold-sdk -DCMAKE_INSTALL_PREFIX=install-dir ..</samp>
3. build and install
<samp>nmake
nmake install</samp>

## Linux
1. clone source code from github
<samp>git clone https://github.com/shihchinw/rlShaders.git source-root</samp>
2. run `cmake` command for configuration:
<samp>cd source-root
mkdir build
cd build
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DARNOLD_ROOT=path-to-arnold-sdk -DCMAKE_INSTALL_PREFIX=install-dir ..</samp>
3. build and install
<samp>make
make install</samp>

# Test Utility

<samp>python testsuite/runtest.py render --ai path-to-MtoA -l rlShaders/shaders --sn 1..10</samp> This command would render test images of test case 1 to 10 and compare the differences against the reference image with `idiff` tool from OpenImageIO. 

Since I only have trial version of Arnold renderer, there are watermarks in all the reference images. Which means the images can't be 100% identical due to the noise in wartermarks. Thus I use RMS error as a metric for comparison instead of the result directly returned from `idiff`.

To view the reference and test images together, we could use `display` sub-command:
<samp>python testsuite/runtest.py display --sn 1 3 5</samp> This will popup an _image viewer_ (`iv` tool from OpenImageIO) loaded with reference/test image, we can use page-up/down to switch the image.

