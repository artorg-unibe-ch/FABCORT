# Bash script to run abaqus simulations and get homogenized stress
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_154_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_161_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_166_M/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_111_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_111_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_112_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_112_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_112_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_112_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_112_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_112_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_113_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_113_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_113_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_113_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_113_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_113_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_114_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_114_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_114_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_114_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_114_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_114_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_121_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_121_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_121_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_121_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_121_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_121_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_122_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_122_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_122_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_122_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_122_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_122_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_123_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_123_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_123_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_123_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_123_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_123_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_124_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_124_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_124_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_124_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_124_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_124_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_211_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_211_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_211_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_211_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_211_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_211_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_212_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_212_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_212_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_212_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_212_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_212_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_213_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_213_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_213_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_213_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_213_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_213_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_214_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_214_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_214_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_214_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_214_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_214_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_221_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_221_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_221_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_221_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_221_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_221_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_222_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_222_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_222_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_222_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_222_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_222_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_223_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_223_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_223_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_223_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_223_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_223_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_224_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_224_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_224_Transverse inp="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_224_Main_Transverse.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_224_Transverse.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_L/ROI_224_Transverse.out"  size="156;156;156" spec="Stress" 
rm * 
abaqus interactive job=ROI_111_Isotropic inp="/home/ms20s284/FABCORT/Homogenization/2012_203_M/ROI_111_Main_Isotropic.inp" cpus=16
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations5/ROI_111_Isotropic.odb"  out="/home/ms20s284/FABCORT/Homogenization/2012_203_M/ROI_111_Isotropic.out"  size="156;156;156" spec="Stress" 
rm * 
