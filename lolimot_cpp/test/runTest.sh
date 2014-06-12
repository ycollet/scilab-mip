#!/bin/bash
echo "Execution du test Cone Noise"
../TrainLolimot2.exe Modele/Modele_Test_Cone_Noise.txt
echo "Execution du test Sinus Noise"
../TrainLolimot2.exe Modele/Modele_Test_Sinus_Noise.txt
echo "Execution du test Temporel Noise"
../TrainLolimot2.exe Modele/Modele_Test_Temporel_Noise.txt
echo "Execution du test Cone woNoise"
../TrainLolimot2.exe Modele/Modele_Test_Cone_woNoise.txt
echo "Execution du test Sinus woNoise"
../TrainLolimot2.exe Modele/Modele_Test_Sinus_woNoise.txt 
echo "Execution du test Temporel woNoise"
../TrainLolimot2.exe Modele/Modele_Test_Temporel_woNoise.txt
echo "Execution du test Plan Noise"
../TrainLolimot2.exe Modele/Modele_Test_Plan_Noise.txt 
echo "Execution du test Temporel Gap Noise"
../TrainLolimot2.exe Modele/Modele_Test_Temporel_Gap_Noise.txt
echo "Execution du test Plan woNoise"
../TrainLolimot2.exe Modele/Modele_Test_Plan_woNoise.txt 
echo "Execution du test Temporel Gap woNoise"
../TrainLolimot2.exe Modele/Modele_Test_Temporel_Gap_woNoise.tx