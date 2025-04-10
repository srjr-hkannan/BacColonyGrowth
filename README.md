# BacColonyGrowth
Code to simulate bacterial colony growth and emergent nutrient gradients  
**Contact:** hkannan@ucsd.edu, paul.sun@csulb.edu

### SOFTWARE DEPENDENCIES:  
+Main agent-based hybrid reaction diffusion model used in  https://doi.org/10.1101/2023.08.27.554977 can be found in the “AgentBasedModel-BacColonyGrowth” folder. 
+Custom code used to analyze microscopy images of colonies can be found in "ImageAnalysis" folder.
+Reduced reaction-diffusion simulation to study dimension dependence of solutions to metabolite concentration can be found in "ReducedReactionDiffusionModel" folder.

### SOFTWARE DEPENDENCIES:  
The agent-based code for colony growth simulations has been tested with the following software versions:  
+	**Operating system**: Ubuntu 20.04 LTS  
+ **Compiler**: g++ (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0 and OpenMP 201511   <br/><br/>
No non-standard hardware is required to run simulations.  

### INSTALLATION GUIDE / INSTRUCTIONS FOR USE:
To compile code:
+ Change current working directory to the “AgentBasedModel-BacColonyGrowth/run_sim” folder.
+ Compilation command : ```g++ ../*.cpp -fopenmp -O3 -o CompiledCode```

Input file for simulation:
+	The input parameters are provided as a “.txt” file (see Supplementary Table of https://doi.org/10.1101/2023.08.27.554977 for description of parameters)
+ Input file for demo (small) simulation: **inputFile_20mM_SmallSimulation.txt**
+ Input file for full simulation of colony dynamics with 20 mM initial glucose concentration (see : https://doi.org/10.1101/2023.08.27.554977):     **inputFile_20mM_LargeSimulation.txt**

To run code:  
+ In the “AgentBasedModel-BacColonyGrowth/run_sim” folder create a new folder named “output” to store output files of simulation.  
+ Execute command: ```./CompiledCode inputFile_20mM_SmallSimulation.txt 4 output/```
  
### DEMO DESCRIPTION
+	The demo simulates *E.coli* colony growth for the first 10 h on a small agar region with 20 mM initial glucose concentration.
+	To run the demo simulation, follow the instructions described in previous sections with the input file **“inputFile_20mM_SmallSimulation.txt”**.
+	The expected output can be found in the “expected_output.zip” folder inside the “run_sim” folder. 
+	The runtime of the demo on a personal laptop with a AMD Ryzen 7 5800H with Radeon Graphics (3.20 GHz) processor and 16 GB RAM was ~1 hour. 

