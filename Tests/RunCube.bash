abaqus interactive job=Cube_Isotropic inp="/home/ms20s284/FABCORT/Test/Cube_Isotropic.inp" cpus=4
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Test/Cube_Isotropic.odb"  out="/home/ms20s284/FABCORT/Test/Cube_Isotropic.out"  size="1;1;1" spec="Stress" 
abaqus interactive job=Cube_Transverse inp="/home/ms20s284/FABCORT/Test/Cube_Transverse.inp" cpus=4
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Test/Cube_Transverse.odb"  out="/home/ms20s284/FABCORT/Test/Cube_Transverse.out"  size="1;1;1" spec="Stress" 
abaqus interactive job=Cube_Fabric inp="/home/ms20s284/FABCORT/Test/Cube_Fabric.inp" cpus=4
abaqus python "/home/ms20s284/FABCORT/Scripts/abqSeReader.py" in="/home/ms20s284/FABCORT/Test/Cube_Fabric.odb"  out="/home/ms20s284/FABCORT/Test/Cube_Fabric.out"  size="1;1;1" spec="Stress" 
