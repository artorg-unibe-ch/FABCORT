# Bash script to run abaqus simulations and get homogenized stress
abaqus interactive job=ROI_111 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_111_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_111.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_111.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_112 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_112_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_112.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_112.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_113 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_113_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_113.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_113.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_114 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_114_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_114.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_114.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_115 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_115_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_115.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_115.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_121 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_121_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_121.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_121.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_122 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_122_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_122.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_122.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_123 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_123_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_123.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_123.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_124 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_124_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_124.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_124.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_125 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_125_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_125.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_125.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_131 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_131_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_131.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_131.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_132 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_132_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_132.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_132.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_133 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_133_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_133.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_133.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_134 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_134_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_134.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_134.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_135 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_135_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_135.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_135.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_211 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_211_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_211.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_211.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_212 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_212_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_212.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_212.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_213 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_213_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_213.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_213.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_214 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_214_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_214.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_214.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_215 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_215_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_215.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_215.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_221 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_221_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_221.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_221.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_222 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_222_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_222.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_222.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_223 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_223_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_223.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_223.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_224 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_224_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_224.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_224.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_225 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_225_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_225.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_225.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_231 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_231_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_231.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_231.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_232 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_232_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_232.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_232.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_233 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_233_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_233.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_233.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_234 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_234_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_234.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_234.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_235 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_235_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_235.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_235.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_311 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_311_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_311.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_311.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_312 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_312_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_312.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_312.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_313 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_313_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_313.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_313.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_314 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_314_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_314.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_314.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_315 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_315_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_315.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_315.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_321 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_321_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_321.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_321.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_322 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_322_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_322.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_322.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_323 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_323_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_323.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_323.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_324 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_324_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_324.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_324.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_325 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_325_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_325.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_325.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_331 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_331_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_331.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_331.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_332 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_332_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_332.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_332.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_333 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_333_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_333.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_333.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_334 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_334_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_334.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_334.out"  size="154;154;154" spec="Stress" 
rm * 
abaqus interactive job=ROI_335 inp="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_335_Main.inp" cpus=24
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Simulations/ROI_335.odb"  out="/home/ms20s284/FABCORT/ConvergenceStudy/ROI_335.out"  size="154;154;154" spec="Stress" 
rm * 
