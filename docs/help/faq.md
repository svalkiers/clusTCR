---
layout: default
title: FAQ
nav_order: 1
parent: Help
---

##  FAQ

> **Is *clusTCR* available for Windows users?**

*clusTCR* relies upon the Faiss library, which is officially only supported for OSX and Linux. Windows users may install a [WSL](https://ubuntu.com/wsl) in order to properly run *clusTCR*.



> **I wanna use the *clustcr-gpu* version, but I don't know if my device has the required hardware.**

*clustcr-gpu* relies on the [Faiss library](https://github.com/facebookresearch/faiss), which on its turn relies on the cudatoolkit.Only specific GPUs with CUDA-compatibility can be used for *clusTCR*'s GPU function. First, you will need to find out which graphics hardware is available in your device. Depending on your OS, you will need to do one of the following things:

- **Linux**

  ```bash
  $ lspci | grep VGA
  ```

  Alternatively, go to *Settings > About > Graphics*.

- **macOS**

  For Mac users, go to *Apple menu > About this Mac > System Report... > Hardware > Graphics/Displays*.

- **Windows**

  Windows users can find graphics card information through *Start > Device Manager > Display adapters*.

Now that you have identified the type of GPU in your device, you need to check its compatibilities. *clustcr-gpu* is only supported for NVIDIA GPUs with CUDA capabilities of minimum compute capability 3.5 ([source](https://github.com/facebookresearch/faiss/wiki/Faiss-on-the-GPU)). To check if your GPU is compatible, visit the [NVIDIA website](https://developer.nvidia.com/cuda-gpus).