---
title: "rlShaders"
description: ""
image: "images/header/index.jpg"
---

# Paper <=> Code

__rlShaders__ is a project for physically-based rendering experiments in Arnold renderer. It primarily focuses on `importance sampling` of recent PBR shading models. It's implemented for learning purposes, you can easily find the corresponding equations between paper and source code.

This also implies that current implementation is not production ready, the performance might be slower than the native shaders in Arnold. So far, it contains following shaders:

* <span class="orange">_rlGgx:_</span> GGX + Oren–Nayar BRDF
* <span class="orange">_rlDisney:_</span> Disney "principled" BRDF
* <span class="orange">_rlSkin:_</span> Normalized diffusion profile of BSSRDF & two-layer GGX specular

<iframe src="https://player.vimeo.com/video/150344036" width="500" height="281" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe>

# Development Journals

* [Implementing GGX BRDF in Arnold with Multiple Importance Sampling](http://shihchinw.github.io/2015/06/implementing-ggx-brdf-in-arnold-with-multiple-importance-sampling.html)
* [Implementing Disney Principled BRDF in Arnold](http://shihchinw.github.io/2015/07/implementing-disney-principled-brdf-in-arnold.html)
* [Sampling Visible Normals for GGX BRDF](http://shihchinw.github.io/2015/08/sampling-visible-normals-for-ggx-brdf.html)
* [BSSRDF Importance Sampling of Normalized Diffusion](http://shihchinw.github.io/2015/10/bssrdf-importance-sampling-of-normalized-diffusion.html)
* [Realistic Human Skin with Normalized Diffusion & GGX](http://shihchinw.github.io/2015/12/realistic-human-skin-with-normalized-diffusion-ggx.html)
