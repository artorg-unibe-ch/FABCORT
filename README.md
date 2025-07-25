# FABCORT: Fabric and Elasticity in Cortical and Trabecular Bone

This repository contains Python scripts to process, analyze, and visualize the mechanical and structural properties of bone microarchitecture, focusing on cortical and trabecular regions. It includes tools for ROI extraction, homogenization, fabric tensor analysis, and micromechanical model fitting (Zysset–Curnier).

---

## 📁 Repository Structure

```
├── Scripts/                # Main processing scripts
├── Data_Cortical/          # Input cortical scan files (.mhd)
├── Data_Trabecular/        # Input trabecular data and .fab files
├── Results/                # Outputs including ROIs, plots, fitted models, etc.
│   ├── Cortical/
│   ├── Tests/
│   └── Fabric.png, CV_Rho.png, *.tex, *.csv
├── Utils/                  # Utility modules (Time, Tensor, Read, Image)
└── README.md               # This file
```
>Note: Folders starting with 99 are temporary during repo restructuration
---

## 🧠 Main Functionalities

### 🖼 1. Scan Visualization
**Script:** `PlotScan.py`  
Plots 3D medical image scans using PyVista and saves them as screenshots.

---

### 🧩 2. ROI Extraction
**Script:** `ComputeROIs.py`  
- Computes Otsu threshold across samples.
- Extracts cubic ROIs (2×2×4 grid).
- Saves Medtool-compatible ROI files.
- Optionally plots ROI layout.

---

### 🛠 3. Abaqus Simulation Setup
**Script:** `GenerateAbaqusInputs.py`  
- Writes `.inp` files with tensile and shear steps.
- Generates bash scripts to batch-run Abaqus simulations.
- Handles both isotropic and transversely isotropic cases.

---

### 📈 4. Homogenization Data Extraction
**Script:** `ExtractStiffness.py`  
- Reads simulation stress data and computes stiffness tensors.
- Rotates tensors into fabric-aligned coordinate system.
- Outputs isotropic and transverse datasets as `.csv`.

---

### 📊 5. Fabric Tensor Projection
**Script:** `PlotFabric.py`  
- Reads MIL-based fabric tensors from `.fab` files.
- Projects ellipsoids onto anatomical planes.
- Plots and compares cortical vs trabecular structure.

---

### 🧮 6. Coefficient of Variation Analysis
**Script:** `ComputeCV.py`  
- Calculates local CV in bone volume fraction across subregions of each ROI.
- Compares trabecular vs cortical variability.
- Outputs a density vs CV plot.

---

### 🧪 7. Zysset–Curnier Model Fitting
**Script:** `FitZyssetCurnier.py`  
- Fits homogenized stiffness to the Zysset–Curnier model (log-log regression).
- Computes confidence intervals, adjusted R², NE.
- Generates plots and LaTeX tables of model parameters.

---

## ⚙️ Requirements

Make sure you have the following Python packages installed:

```bash
numpy
pandas
matplotlib
SimpleITK
pyvista
scipy
```

To install them, run:

```bash
conda env create -f Requirements.yml
```

>Note: Requirement.yml is not yet created

---

## 📌 Notes & Next Steps

- ✅ Check trabecular bone Zysset–Curnier regression to match results from the original publication.
- 🧾 Finish documenting `YandandCowinModel.py`, including:
  - A proper module docstring
  - Export of fitted parameters to a LaTeX table
- 🔬 Implement two modified Zysset–Curnier models for transverse isotropy (see `FabCortMathematica.pdf`):
  - Model with 3 constants and 4 exponents
  - Model with 3 constants and 7 exponents
- 🔬 Implement comparison between RUS measurement and numerical simulation of cortical bone with isotropic and transverse isotropic matrix
---

## 📚 References

1. Zysset, P.K. & Curnier, A. (1995).
*An alternative model for anisotropic elasticity based on fabric tensors*.
Mechanics of Materials, 21(4), 243–250.
https://doi.org/10.1016/0167-6636(95)00018-6

2. Panyasantisuk, J. et al. (2015).
*Comparison of Mixed and Kinematic Uniform Boundary Conditions in Homogenized Elasticity of Femoral Trabecular Bone Using Microfinite Element Analyses*.
Journal of Biomechanical Engineering, 137(1).
https://doi.org/10.1115/1.4028968

3. Yang, G. et al. (1999).
*The Anisotropic Hooke's Law for Cancellous Bone and Wood*.
Journal of Elasticity, 53, 125–146.

4. Mehrabadi, M.M. & Cowin, S.C. (1990).
*Eigentensors of Linear Anisotropic Elastic Materials*.
QJMAM.
http://qjmam.oxfordjournals.org/

5. Cowin, S.C. & Yang, G. (1997).
*Averaging Anisotropic Elastic Constant Data*.
Journal of Elasticity, 46, 151–180.

6. Cai, X. et al. (2019).
*Anisotropic elastic properties of human femoral cortical bone and relationships with composition and microstructure in elderly*.
Acta Biomaterialia, 90, 254–266.
https://doi.org/10.1016/j.actbio.2019.03.043

7. Simon, M. et al. (2022).
*Fabric-elasticity relationships of tibial trabecular bone are similar in osteogenesis imperfecta and healthy individuals*.
Bone, 155, 116282.
https://doi.org/10.1016/j.bone.2021.116282 

8. Franzoso, G. & Zysset, P. (2009).
*Elastic anisotropy of human cortical bone secondary osteons measured by nanoindentation*.
J Biomech Eng, 131(2).
https://api.semanticscholar.org/CorpusID:25765365

9. Enrico, D., Schmidt, R. & Zysset, P. (2012).
*Microindentation can discriminate between damaged and intact human bone tissue*.
Bone, 50(4).
https://api.semanticscholar.org/CorpusID:23349859  

