# Ptycho_prep
Scripts to prep data and analyze output of transmission Ptycho measurements

Requirements: Compiled Ptycholib and Matlab

Procedure:

Step 1: Fire up Matlab, generate a probe guess in Matlab. 

	      probe=v2_quick_probe(0.1,0,0,0,10.4,800000,1)
	      
          --parameters are thickness, theta, twotheta, defocus, E (kev), detector distance, 1)
	  
        Save it, use probe_test.ipynb to look at the probe and save to text file

Step 2: Use the text file from probe_test.ipynb as input probe guess

        ---Alternatively, use the probe guesses in 2018R1/20180402/ but this is on 512 array!

Step 3: Run setup_ptycho to generate h5 files from tif files. 

	      Remember to change imgsperpt and imstart if exits.

Step 4: Modify runptycho.sh appropriately and run

