abaqus interactive job=ROI_1_F1 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_1_F1.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F1.out"  size="155;155;155" spec="Stress" 
rm * 
abaqus interactive job=ROI_1_F2 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_1_F2.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F2.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_1_F4 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_1_F4.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_1_F4.out"  size="158;158;158" spec="Stress" 
rm * 
abaqus interactive job=ROI_2_F1 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_2_F1.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F1.out"  size="155;155;155" spec="Stress" 
rm * 
abaqus interactive job=ROI_2_F2 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_2_F2.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F2.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_2_F4 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_2_F4.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_2_F4.out"  size="158;158;158" spec="Stress" 
rm * 
abaqus interactive job=ROI_3_F1 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F1_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_3_F1.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F1.out"  size="155;155;155" spec="Stress" 
rm * 
abaqus interactive job=ROI_3_F2 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F2_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_3_F2.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F2.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_3_F4 inp="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F4_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_3_F4.odb"  out="/home/ms20s284/FABCORT/ResolutionEffect/ROI_3_F4.out"  size="158;158;158" spec="Stress" 
rm * 
