abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_111_Main_Transverse.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
