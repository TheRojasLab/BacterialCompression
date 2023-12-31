# BacterialCompression
https://www.biorxiv.org/content/10.1101/2022.08.12.503793v1

# Mechanical Compression Induces Persistent Bacterial Growth 

## Description

The scripts ImageAlign, BTFluo, BacTrack2 are for analyzing fluorescent _E. coli_ cells in a CellASIC perfusion chip using both phase and fluorescence channels. They will measure variables
such as fluorescence and growth rate.
The above scripts with _ROI function in the same manner but are for analyzing windows of the whole perfusion chamber to be combined with the script Combining_uM_Spaces.
The script ZStackHeight will generate a graph of the height of the CellASIC perfusion chamber. This requires a tiled image of the perfusion chamber perfused with a fluorescent dye to generate
the height gradient. 



### Executing program

* How to run the program
* Have three directories containing three image stacks, one with the GFP channel, one with the phase channel and one containg a cropped sequence of images without any cells for the script to use as a
* reference to correct for drift.
* Enter the directory locations in ImageAlign and run the script. You will create two new directories 1_a_ and 2_a_
* Enter directory locations in BacTrack2 and run the script. Two files will be created _BT and _BTlab
* Enter location of the _BT and _BTlab files in BTfluo and run the script. A file will be created called _BTfluo that contains your data

* To run ZStackHeight requires an image of the entire stitched perfusion chamber containg fluorescence. Enter trappos when _E. coli_ cells begin to be trapped, this is in uM on a scale from 0 -> 1500
* e.g. 950. This sets the value of the height at this spot at 1uM.
  
* The ROI scripts require identification of analysisregion. Specify either 200/400/600/800/1000/1200/1400. Once you have all of these in BTFluo format. Use Combining_uM_spaces to collate. 
```


## Authors

Contributors names and contact info

Guy Mason (guymason@nyu.edu)
Enrique Rojas (rojas@nyu.edu)


