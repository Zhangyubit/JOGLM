## Overview

This repository provides MATLAB implementations of **joint time–vertex graph filter banks**, including **oversampled graph filter banks (OGFB)**, for applications such as graph signal denoising and reconstruction.

The proposed framework extends classical graph filter bank designs to the **joint time–vertex domain**, enabling multiresolution analysis for signals defined on product graphs.

---

## Dependencies

This code **must be executed within the GSPBox framework**:

* GSPBox

Please install and initialize GSPBox before running the code:

```matlab
gsp_start;
```

---

## References

Parts of this implementation are inspired by the following foundational works on graph filter banks:

* S. K. Narang and A. Ortega, Perfect Reconstruction Two-Channel Wavelet Filter Banks for Graph Structured Data,
  *IEEE Transactions on Signal Processing*, 2012.

* Y. Tanaka and A. Sakiyama, M-Channel Oversampled Graph Filter Banks, *IEEE Transactions on Signal Processing*, 2014.

* A. Sakiyama and Y. Tanaka, Oversampled Graph Laplacian Matrix for Graph Filter Banks, *IEEE Transactions on Signal Processing*, 2014.

Readers are encouraged to consult these papers for theoretical background if needed.

---

## Copyright

© 2025–2026 **Yu Zhang**. All rights reserved.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{zhang2026twochannelfilterbanksjoint,
  title={Two-Channel Filter Banks on Joint Time-Vertex Graphs with Oversampled Graph Laplacian Matrix}, 
  author={Yu Zhang and Bing-Zhao Li},
  year={2026},
  eprint={2511.11768},
  archivePrefix={arXiv},
  primaryClass={math.GM},
  url={https://arxiv.org/abs/2511.11768}
}
```

---

## Notes

* The implementation focuses on:

  * Joint time–vertex filter bank design
  * Oversampled graph Laplacian construction
  * Denoising and reconstruction experiments

* Example scripts are provided for:

  * Synthetic graph signals
  * Image data
  * Video sequences
