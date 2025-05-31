# Colony-Exit-Bat-Simulation

BatModel– Readme (May 2025)
This guide explains how to run and modify the bat navigation simulation described in the article using the accompanying MATLAB GUI. 
1. Installation Instructions
1.1. Extract the contents of the .rar archive to a local folder (e.g., C:\BatModel\).

2. Launch the MATLAB GUI
2.1. Open MATLAB®.
2.2. In the Command Window, set the current path to the folder where the files were extracted: (e.g., cd(‘c:\BatModel\’)
2.3. Run the GUI by typing: BatGUI1
2.4. A Graphical User Interface that controls and executes the simulation will be opened:



 


3. Running the Simulation (1)
3.1. Load predefined parameters from an Excel file (DefaultParamsTable.xlsx in the \DATA\ folder) or define your own .xlsx file following the same structure.
3.2. You can also adjust key simulation parameters directly in the GUI (see below).
3.3. To start the simulation, press Run.
3.4. Press Send to WS to export all simulation data to the MATLAB workspace for additional analysis (see below).

4. GUI Parameter Settings 
4.1. Main Parameters (2)
4.1.1. Set the number of bats (1–20) and prey items (use zero for cave exit simulations).
4.1.2. Define the total simulation time (default = 15 seconds).
4.1.3. Leave other settings unchanged unless specifically needed.

4.2. Bat Behavior (3)
4.2.1. Confirm the Mode is set to caveExit (default).
4.2.2. Ensure Avoid Obstacles is set to on.
4.2.3. Set memory size for multi-call integration via Obs. Mem (default = 5 calls).
4.2.4. You may activate the confusion mode and toggle multi-call clustering, as discussed in the article.

4.3. Bat Flight (4)
4.3.1. Adjust flight speed and maximum acceleration to control maneuverability.
4.3.2. Other flight parameters can remain at default values.

4.4. Bat Sonar (Echolocation) (5)
4.4.1. Basic sonar parameters can be modified via the GUI.
4.4.2. For advanced echolocation behavior, edit the parameter Excel sheet.

5. Plotting Results
5.1. Basic Visualization (6,7)
5.1.1. Select individual or all bats to visualize.
5.1.2. Choose whether to plot within the GUI or in a new MATLAB figure (new Fig checkbox).
5.1.3. Use Plot X-Y (upper axis) and Plot time-Graph (LOWER AXIS) to explore spatial and temporal data.

5.2. Advanced Analysis (8)
5.2.1. Acoustics Rx: Visualizes received signals of the selected bat over time (signal type by color). 
5.2.2. Analyze Call: 
Analyze Rx / Analyze Man: Examine specific calls, their detections, and the bat’s behavioral decisions during those calls.

6. Output Data for Analysis
6.1. Simulation results are stored in the BatDATA structure. Use Send to WS to export data to the MATLAB workspace or save/load .mat files.
6.2. BatDATA consists the following fields:
BatDATA = struct with fields:
                       	PREY: [1×(number_of_prey_items) struct]
                        	BAT: [1×( number_of_bats) struct]
                  	AllParams: [1×1 struct]
  		FlightInterferceSummary: [1×1 struct]
              		FilterBank: [1×1 struct]

6.3. Main fields in BatDATA to analyze:
6.3.1. BatDATA.AllParams – All parameter values for the simulation
6.3.2. BatDATA.FlightInterferceSummary – Summary of masking and jamming events.
6.3.3. BatDATA. BAT(x). InterReportStrctOnLine – Per-bat trajectory and detection details.

7. Core Files and Examples
7.1.1. BatGUI1.m – Main GUI launcher.
7.1.2. BatFlightForGui.m – Core simulation function.
7.1.3. \DATA\ DefaultParamsTable_CaveExit_FInal_PK.xlsx – Default parameter file.
7.1.4. \BatDATA_output\ BatData_.mat – Example of simulation output 
7.1.5. \Experiments Code\*.m - Batch-run examples and automated testing scripts.


Contact
Omer Mazar
Omer_mazar@yahoo.com
Updated: May 2025


