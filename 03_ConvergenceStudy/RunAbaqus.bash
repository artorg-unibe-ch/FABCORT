# Bash script to run abaqus simulations and get homogenized stress
abaqus interactive job=ROI_S10_N10 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N10_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N10.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N10.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N6 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N6_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N6.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N6.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N7 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N7_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N7.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N7.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N8 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N8_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N8.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N8.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S10_N9 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N9_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S10_N9.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S10_N9.out"  size="153;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S1_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S1_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S1_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S1_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S2_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S2_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S2_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S2_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S2_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S2_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S2_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S2_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S3_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S3_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S3_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S3_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S3_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S3_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S3_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S4_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S4_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S4_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S4_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S4_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S4_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S4_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S4_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S4_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S5_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S5_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S5_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S5_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S5_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S5_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S5_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S5_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S5_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S5_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S5_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S6_N6 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N6_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S6_N6.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S6_N6.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N6 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N6_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N6.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N6.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S7_N7 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N7_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S7_N7.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S7_N7.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N6 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N6_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N6.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N6.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N7 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N7_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N7.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N7.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S8_N8 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N8_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S8_N8.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S8_N8.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N1 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N1.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N1.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N2 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N2.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N2.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N3 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N3_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N3.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N3.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N4 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N4.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N4.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N5 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N5_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N5.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N5.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N6 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N6_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N6.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N6.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N7 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N7_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N7.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N7.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N8 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N8_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N8.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N8.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_S9_N9 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N9_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_S9_N9.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_S9_N9.out"  size="154;154;154" spec="Stress" 
rm * 
