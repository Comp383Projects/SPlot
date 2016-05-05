# SPlot 2.0:
SPlot is a tool used to find HGT events and RUCPs throughout genomes built in JAVA with the Glimpse data visualization framework and JeUtils BLAST searching package (all packaged in SPlot JAR).

#System Requirements:
Java Version 8 <br>
OpenGL Compatible Graphics Card/Drivers

#How to Download:
Download is availible via Google Drive: https://drive.google.com/file/d/0B-yWo52TprPralQxbU5JNlRlQnM/view?usp=sharing

#How to Setup SPlot run:
1) Double click on SPlot JAR file. <br>
2) Type in kmer size (Note: large kmer sizes will not produce informative SPlots). Recommended: 3 or 4. <br>
3) Type in window size (Note: the limit is system dependent -> size 500-5000 recommended, with a larger genome requiring larger size). <br>
4) Type in Sequence Names (these are used for output file naming, as well as labeling the SPlot) <br>
5) Click on the First/Second Sequence button and navigate to your sequence (NOTE: FASTA files do not work yet-> only sequences allowed, aka: all headers must be deleted). <br>
6) Submit for SPlot Visualization. (Note: currently no loading screen) <br>

#Deriving Results from SPlot
With use of a mouse, you are able to use the scroll wheel to zoom in and out of SPlot and look at the correlations determined by Pearsons' Correlation Coefficient. <br>
In the backend GUI, several options exist: <br>
1) Individual Zones: Select either the X sequence or Y sequence, and enter the Start/Stop Coordinates of that sequence. From there, you can a) Write the Sequence to file or b) BLAST the sequence against NCBI's Sequence Database. All files will be saved in SPlot's working directory. <br>
2) Whole SPlot Analysis: This will identify rows or columns with a <.5 average Pearson Correlation Coefficient and <.7 maximum Pearson Correlation Coefficient. Then, these rows/columns will be a) written out to file in FASTA format (FileAnalysis) or b) BLASTed using JeUtils and written to file (note, this should only be run in BLAST non-peak hours-> nights and weekends).
