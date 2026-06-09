# Custom analysis code — Dynamic remodeling of the pancreas immune landscape in obesity

Alexey Koshkin¹, Kranthi Kiran Kishore Tanagala¹, Anna Eichinger²·³, Michael Chait¹, Aoife Young¹, Shanila Shakil¹, Junichi Yoshikawa⁴, Yosuke Sakamoto⁴, Steven B. Wells⁴, Ladan Fazlollahi¹, Xiaojuan Chen⁵, Boris Reizis², Donna L. Farber⁴·⁶ and Stuart P. Weisberg¹·*

1. Department of Pathology and Cell Biology, Columbia University Irving Medical Center, New York, NY 10032

2. Department of Pathology, NYU Grossman School of Medicine, New York, NY 10016

3. Department of Pediatrics, Dr. von Hauner Children's Hospital, University Hospital, Ludwig-Maximilians-University (LMU), Munich, Germany

4. Department of Microbiology and Immunology, Columbia University Irving Medical Center, New York, NY 10032

5. Columbia Center for Translational Immunology, Columbia University Irving Medical Center, New York, NY 10032

6. Department of Surgery, Columbia University Irving Medical Center, New York, NY 10032

*Corresponding author: Stuart P. Weisberg

Custom scripts written for this study. These are the analyses developed specifically for this work; all other analyses used published tools as cited in the Methods. This repository contains only the custom code.

## Contents

Spatial neighborhood analysis (Figure 2):

- Fig2_c_d_e_S2B.ipynb — generates Figure 2c, 2d, 2e and Figure S2B.

- Fig2_GH_CD11c_hi_lo_figure.ipynb — generates Figure 2G and 2H (CD11c high/low).

- vectra_lib.py, vectra_lib_PCF.py — helper functions imported by the Figure 2 notebooks (multiplex imaging / neighborhood analysis).

Supplementary figure scripts:

- S3D_S3E_adipose_Tcell.R — pancreatic T/NK gene programs scored on the Massier adipose lymphoid atlas (GSE225700); panels S3D, S3E.

- S5A_S5B_adipose_myeloid.R — pancreatic macrophage gene programs scored on the Massier adipose myeloid atlas (GSE225700); panels S5A, S5B.

- S8_CD163_histogram.R — CD163 intensity across macrophage subsets vs. acinar/islet controls, with a 5th-percentile threshold; panel S8 (CD163).

- S8_GZMB_density.R — GZMB intensity density for macrophages vs. T cells, GZMB+ defined as raw intensity > 2; panel S8 (GZMB).

## Data

The adipose scripts use the Massier et al. adipose atlas, publicly available at GEO GSE225700; the input objects are the lymphoid and myeloid subsets of that data.

The Figure 2 and S8 scripts use imaging data from this study. The underlying objects contain donor-derived data and are not included here; they are available from the corresponding author on reasonable request. Per-cell values for the figure panels are provided in the paper's Source Data.

Set the input/output paths at the top of each script before running.

## Code vs. published panels

- S3D / S3E — the script scores the full gene set including a 25-gene obesity-IFN program; the published panels show the three-program subset (cytotoxic function, chemotaxis, negative regulation). The CCL5/TCF7/TNFAIP3 UMAPs on the S3D page were generated separately.

- S5A — matches the published panel as written.

- S5B — the script produces a 10-program × 16-cluster heatmap; the published panel collapses the 16 myeloid clusters into three groups (fetal-like, LAM, monocytic), assembled manually.
