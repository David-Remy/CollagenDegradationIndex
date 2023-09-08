/*  Tested on ImageJ 1.53c / Windows
 *  Created in March 2021 by David Remy (Main Author) and Anne-Sophie Macé (Advices & corrections)

 *  This macro computes a degradation index, which is the ratio between the area of collagen cleavage and the number of nuclei in the image
 *  Batch quantification included: input should be a folder containing all acquisitions
 *  
 *  Are expected (at least) those 2 stainings: collagen 1-¾C & nucleus staining
 */

 // Clears everything
roiManager("reset");
run("Clear Results");
run("Close All");
setOption("ExpandableArrays", true); // In ImageJ 1.53g and later, arrays automatically expand in size as needed. This option is used in case of early versions of ImageJ
if( isOpen("Summary") ){
	selectWindow("Summary");
	run("Close");
}
print("\\Clear");

/////////////////////////////////////////////////////////////
////// begining of parameters customozible by the user //////
/////////////////////////////////////////////////////////////
// size of the rolling ball in pixels for the substract bg of the collagen 3/4 
rolling_ball_param_col = 5;
// minimum area in pixels for Analyze particle in the collagen 3/4 detection
analyze_particle_min_col = 10;
// minimum area in pixels of the nuclei to detect
nuclei_min_area = 1000;
// extensions of the file to study [tested on ND only]
ext_file = ".tif";
/////////////////////////////////////////////////////////////
//////// end of parameters customozible by the user /////////
/////////////////////////////////////////////////////////////
 

// Select input directory and create an Analysis subfolder
dir_input=getDirectory("Select input directory");
dir_output=dir_input+"Collagen_Degradation_Analysis"+File.separator; 
if( !File.exists(dir_input+"Collagen_Degradation_Analysis") ) // creates directory if does not exist
	File.makeDirectory(dir_output);

// Dialog box for channel representation and the choices of slices for z-projection
Dialog.create("Channel representation");
Dialog.addNumber("Nucleus channel", 1);
Dialog.addNumber("Col1-3/4C channel", 3);
Dialog.addCheckbox("If z-stack, Max Projection on the entire z-stack ?" , false);
Dialog.show();
Nuc_channel=Dialog.getNumber();
Degra_channel=Dialog.getNumber();
proj=Dialog.getCheckbox();

//list every files names within the directory
Filelist=getFileList(dir_input);
Array.sort(Filelist);

//Initializes arrays to compile all names of images, number of nuclei per image, area of degradation per image and area/numb of nuclei
Name=newArray();
NumNuclei=newArray();
AreaDegradation=newArray();
normArea=newArray();

// loop on all files of the folder
file_treated = 0;
for (i_file=0; i_file<lengthOf(Filelist); i_file++) {
	if(indexOf(Filelist[i_file], ext_file)>0){ // treats each ND file
		shortTitle = substring(Filelist[i_file],0,lastIndexOf(Filelist[i_file],"."));
		
		run("Bio-Formats Importer", "open=["+dir_input+Filelist[i_file]+"] autoscale color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		tit_img = getTitle();
		
		// convert scale in pixels
		run("Set Scale...", "distance=1");
		getDimensions(width, height, channels, slices, frames);
		
		// check that the specified channels exist
		if( channels < maxOf(Degra_channel,Nuc_channel) )
			print("Image "+Filelist[i_file]+" not treated: at least one of the specified channel do not exist");
		
		else{
			print("Treating image "+Filelist[i_file]);
			
			if( !File.exists(dir_output+"Results_"+shortTitle+".xls") ){ // the image was not treated
				// duplicates the channels of interest
				selectWindow(tit_img);
				run("Duplicate...", "title=col1-3/4C duplicate channels="+Degra_channel);
				selectWindow(tit_img);
				run("Duplicate...", "title=nucleus duplicate channels="+Nuc_channel);
				
				//counts the number of nuclei in the field on the z-max projection
				selectWindow("nucleus");	
				if (nSlices >1 ) { // If the image is a stack: maximum z-projection
					run("Z Project...", "projection=[Max Intensity]");
					close("nucleus");
					selectWindow("MAX_nucleus");
					rename("nucleus");
				}
				run("Gaussian Blur...", "sigma=2");
				setAutoThreshold("Default dark");
				run("Analyze Particles...", "size="+nuclei_min_area+"-Infinity exclude display"); //The "exclude" command excludes nuclei that are on the borders of the image
				
				// Nuclei = number of found particules; a dialog box asks the user to confirm or modify this number
				nuclei=nResults;
				NumNuclei[file_treated]=getNumber("Nuclei detected", nuclei);
				
				selectWindow("col1-3/4C");	
				first_plan = 1;
				last_plan = nSlices;
				
				// image is a stack and the user asked to choose the plans of analysis
				if( nSlices > 1 && !proj) {
					// asks the user the z plans on which the projection should be performed
					Dialog.createNonBlocking("Plans choice");
					Dialog.addMessage("Choose the plan to analyze, between 1 and "+slices+" (put the same value if you want only one image)."); 
					Dialog.addMessage("If several plans, a projection will be performed.");
					Dialog.addNumber("First plan (above 1)", first_plan);
					Dialog.addNumber("Last Plan (below "+slices+")", last_plan);
					Dialog.show();
				
					first_plan = Dialog.getNumber();
					last_plan = Dialog.getNumber();
					
					if( first_plan < 1 || first_plan > slices || last_plan < 1 || last_plan > slices)
						print("You choose an uncorrect number for first or last plan, the maximum projection is computed on all slices");
			
					// apply maximum z-projection on the slices specified by the user
					selectWindow("col1-3/4C");
					run("Z Project...", "start="+first_plan+" stop="+last_plan+" projection=[Max Intensity]");
					close("col1-3/4C");
				}
				// Max Projection on all the stack if the user chose that option 
				if( nSlices > 1 && proj) 
					run("Z Project...", "start="+first_plan+" stop="+last_plan+" projection=[Max Intensity]");	
				else // to have the same name
					rename("MAX_col1-3/4C");
				
				selectWindow("MAX_col1-3/4C");
				// pre-processing - user defined threshold - analyze particle
				run("Subtract Background...", "rolling="+rolling_ball_param_col);
				run("Threshold...");
				setAutoThreshold("Default dark");
				waitForUser("Please select the appropriate thresold. Then, press OK !");
				run("Analyze Particles...", "size="+analyze_particle_min_col+"-Infinity clear summarize ");
				
				// saves updated summaries and thresholds 
				close("*");
				selectWindow("Summary");
				IJ.renameResults("Summary","Results");
				Name[file_treated]=shortTitle;
				AreaDegradation[file_treated]=getResult("Total Area", 0);
				normArea[file_treated]=AreaDegradation[file_treated]/NumNuclei[file_treated];
				setResult("NucleiNumber",0,NumNuclei[file_treated]); // to be able to read again results if the macro crashes/the user had to stop
				saveAs("Results", dir_output+"Results_"+shortTitle+".xls");
				file_treated++;
				close("Results");
			}
			else{ // the image was aready treated: we load the results
				print("Results file existed: loaded");
				run("Results... ","open=["+dir_output+"Results_"+shortTitle+".xls]");
				Name[file_treated]=shortTitle;
				AreaDegradation[file_treated] = getResult("Total Area", 0);
				NumNuclei[file_treated] = getResult("NucleiNumber", 0);
				normArea[file_treated] = AreaDegradation[file_treated]/NumNuclei[file_treated];
				file_treated++;
				close("Results");
			}
		}
	}	
}

// final result table for all images of the folder is created and saved
run("Clear Results");
for (i_results = 0; i_results < lengthOf(NumNuclei); i_results++) {
	setResult("Name", i_results, Name[i_results]);
	setResult("Nuclei", i_results, NumNuclei[i_results]);
	setResult("Total Area", i_results, AreaDegradation[i_results]);
	setResult("Area/Nuclei", i_results, normArea[i_results]);
}	
saveAs("Results", dir_output+File.getName(dir_input)+"_NormalizedArea.xls");

run("Close All");

