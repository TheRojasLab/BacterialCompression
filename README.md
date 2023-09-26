# BacterialCompression
https://www.biorxiv.org/content/10.1101/2022.08.12.503793v1

# Mechanical Compression Induces Persistent Bacterial Growth 

## Description

The scripts ImageAlign, BTFluo, BacTrack2 are for analyzing fluorescent _E. coli_ cells in a CellASIC perfusion chip using both phase and fluorescence channels. They will measure variables
such as fluorescence and growth rate.


### Executing program

* How to run the program
* Have three directories containing three image stacks, one with the GFP channel, one with the phase channel and one containg a cropped sequence of images without any cells for the script to use as a
* reference to correct for drift.
* Enter the directory locations in ImageAlign and run the script. You will create two new directories 1_a_ and 2_a_
* Enter directory locations in BacTrack2 and run the script. Two files will be created _BT and _BTlab
* Enter location of the _BT and _BTlab files in BTfluo and run the script. A file will be created called _BTfluo that contains your data

* 
```


## Authors

Contributors names and contact info

Guy Mason (guymason@nyu.edu)
Enrique Rojas (rojas@nyu.edu)


