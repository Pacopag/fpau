fpau
====

Simple model of global land use dynamics

Usage
====

Unzip the files to some folder on your hard drive.
cd to that folder.

Compile (e.g. FPAU.java) with
javac -cp "./src:./lib/jcommon-1.0.17.jar:./lib/jfreechart-1.0.14.jar" -d ./bin ./src/FPAU.java 

There are several ways to run it

1. Run the baseline scenario
java -cp "./bin:./lib/jcommon-1.0.17.jar:./lib/jfreechart-1.0.14.jar" FPAU

2. Run the "splice" to freely choose carrying yield and consumption capacities (last two args are Ky and Kc)
java -cp "./bin:./lib/jcommon-1.0.17.jar:./lib/jfreechart-1.0.14.jar" FPAU splice 3600 2550

3. Run the phase diagram calculation
java -cp "./bin:./lib/jcommon-1.0.17.jar:./lib/jfreechart-1.0.14.jar" FPAU phases

4. Run the sensitivity analysis (last two args are percentage to vary from baseline values and number of steps)
java -cp "./bin:./lib/jcommon-1.0.17.jar:./lib/jfreechart-1.0.14.jar" FPAU sensitivity 10 100

In all cases any data files generated will be created in the data/ folder
