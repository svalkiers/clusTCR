## clusTCR: a Python interface for rapid clustering of large sets of CDR3 sequences

clusTCR is a two-step clustering approach that combines the speed of the [faiss](https://github.com/facebookresearch/faiss) library with the accuracy of [Markov Clustering Algorithm](https://micans.org/mcl/). Compared to other state-of-the-art clustering algorithms ([GLIPH2](http://50.255.35.37:8080/),  [iSMART](https://github.com/s175573/iSMART) and [tcrdist](https://github.com/kmayerb/tcrdist3)), clusTCR show comparable clustering quality, but provides a steep increase in speed and scalability. Using a standard laptop, clusTCR can cluster 1 million CDR3 sequences in under 5 minutes (Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz, using 8 CPUs).

<p align="center">
  <img src="img/workflow.png" alt="drawing" width="800"/>
</p>

