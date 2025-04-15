# Bash script to run abaqus simulations and get homogenized stress
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_067_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_121_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_206_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_211_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_234_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations3/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2010_239_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
