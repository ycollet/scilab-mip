#!/bin/bash
echo "Execution du test Cone Noise"
../lolimot Modele/Modele_Test_Cone_Noise.txt
echo "Execution du test Sinus Noise"
../lolimot Modele/Modele_Test_Sinus_Noise.txt
echo "Execution du test Temporel Noise"
../lolimot Modele/Modele_Test_Temporel_Noise.txt
echo "Execution du test Cone woNoise"
../lolimot Modele/Modele_Test_Cone_woNoise.txt
echo "Execution du test Sinus woNoise"
../lolimot Modele/Modele_Test_Sinus_woNoise.txt 
echo "Execution du test Temporel woNoise"
../lolimot Modele/Modele_Test_Temporel_woNoise.txt
echo "Execution du test Plan Noise"
../lolimot Modele/Modele_Test_Plan_Noise.txt 
echo "Execution du test Temporel Gap Noise"
../lolimot Modele/Modele_Test_Temporel_Gap_Noise.txt
echo "Execution du test Plan woNoise"
../lolimot Modele/Modele_Test_Plan_woNoise.txt 
echo "Execution du test Temporel Gap woNoise"
../lolimot Modele/Modele_Test_Temporel_Gap_woNoise.txt

