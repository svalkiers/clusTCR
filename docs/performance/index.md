---
layout: default
title: Performance
nav_order: 2
---

## Why should I use clusTCR?

*clusTCR* offers a 50x speed improvement (at 1e6 sequences) over GLIPH2, the second fastest TCR/CDR3 clustering available. At the same time, the clustering quality of *clusTCR* remains comparable to other approaches. Benchmarking revealed that *clusTCR* can cluster one million CDR3 sequences in less than 5 minutes.

<p align="center" style="margin-top: 10px">
  <img src="speed.png" alt="drawing" style="width: 70%"/>
</p>
*clusTCR* implements multi-CPU processing to parallelize and speed up the clustering process. Our algorithm also offers GPU support for the calculation of support. This may additionally speed up the overall runtime, although we did not observe significant speed improvements*.

<sub>* : NVIDIA Quadro T2000 Mobile / Max-Q</sub>