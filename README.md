# SeqFISH pipeline

### Description of the pipeline

Our pipeline is a **Matlab®** and **R** based pipeline that performs all the basic steps of seqFISH and smFISH data analysis, namely :
- Image cleaning and pre-processing.
- Automated spot detection/localisation.
- Automated cell detection and segmentation.
- Image stitching and alignment across multiple rounds.

The whole pipeline can be run in an automated manner if the data quality is high enough.

### Mathematical tools used

Our pipeline relies on several powerful mathematical tools coming from different sub-fields :
- Spot detection is performed using the H-dome based spot detector. We developped a modified version of the Mean-Shift algorithm based on KD trees to increase the computational speed of the pixel clustering step (see manuscript).
- Cell segmentation is performed using the RNA spots : a diffusion based method (NJW clustering method) allows to quickly cluster spots without any prior on the shape and size of the cells.
- Image stitching and alignment is performed using the phase correlation method and Fast Fourier Transform (FFT).
