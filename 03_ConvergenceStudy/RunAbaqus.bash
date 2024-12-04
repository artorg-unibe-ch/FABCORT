# Bash script to run abaqus simulations and get homogenized stress
abaqus interactive job=ROI_017 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_017_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_017.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_017.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_018 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_018_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_018.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_018.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_019 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_019_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_019.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_019.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_020 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_020_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_020.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_020.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_021 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_021_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_021.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_021.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_022 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_022_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_022.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_022.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_023 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_023_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_023.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_023.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_024 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_024_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_024.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_024.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_025 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_025_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_025.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_025.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_026 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_026_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_026.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_026.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_027 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_027_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_027.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_027.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_028 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_028_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_028.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_028.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_029 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_029_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_029.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_029.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_030 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_030_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_030.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_030.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_031 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_031_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_031.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_031.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_032 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_032_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_032.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_032.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_038 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_038_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_038.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_038.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_039 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_039_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_039.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_039.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_040 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_040_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_040.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_040.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_041 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_041_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_041.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_041.out"  size="154;154;154" spec="Stress" 
rm * 
