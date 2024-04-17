/* 
DrosophilaGP Macro
This macro will take the segmented images of germplasm (performed in ilastik) and convert the image stacks into .csv files.
This conversion into .csv files will also be performed after several different versions of 3D Objects Counter plugin being performed on the segmented images.
Make sure that the chosen options for 3D Objects Counter are the same as in the image "3D_Objects_Counter_Settings" in the ImageJ_macros folder on Github.
This process can take some time to finish running on all the images in a folder. 
Each image will have its own folder with subfolders that were just created in the ImageJ macro (and will be stored in the output directory).
*/

//Storing the path to the input files (folder with the HDF5 segmentations) and the output directory (where the analysis will be stored)
inputDir =getDirectory("Choose an input directory"); 
outputDir =getDirectory("Choose an output directory");
setBatchMode(true);

//Reading the files in the input directory 
list = getFileList(inputDir);


// LOOP --> Getting the results for the thresholded image
for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
    run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath1 = outputDir + outputFilename;
	File.makeDirectory(thresholdingOutputPath1);
	thresholdingOutputPath2 = outputDir + outputFilename + "/ZSlices/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);
	
	slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);
	
	run("Image Sequence... ", "dir=" + thresholdingOutputPath2 + " format=TIFF digits=2");
	resultsOutputPath = outputDir+outputFilename+"/ZResults/";
	File.makeDirectory(resultsOutputPath);
    zSliceList = getFileList(thresholdingOutputPath2);
    for(j=0;j<zSliceList.length;j++){
        open(thresholdingOutputPath2+zSliceList[j]);
        run("Image to Results"); //Results files are saved as Results0, Results1, ..
        resultsFileName = "Results"+j;
        saveAs("Results",resultsOutputPath+resultsFileName + ".csv");
        close();
    }
    close("*");
    selectWindow("Results"); 
    run("Close");
}

// LOOP --> Getting the results for the image with edges EXCLUDED
for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
    run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath1 = outputDir + outputFilename;
	File.makeDirectory(thresholdingOutputPath1);
	thresholdingOutputPath2 = outputDir + outputFilename + "/CentroidData/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);
	
	slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);
	
	run("3D Objects Counter", "threshold=128 slice=11 min.=4 max.=99999999999 exclude_objects_on_edges statistics");
	CentroidDatafilename = "3D_ObjectCount_ResultsTable.csv";
	saveAs("Results", thresholdingOutputPath2 + CentroidDatafilename);
    close("*");

}

// LOOP --> Getting the results for the image with edges INCLUDED
for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
	run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath2 = outputDir + outputFilename + "/CentroidData_All/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);
	
    slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);
	
	run("3D Objects Counter", "threshold=128 slice=11 min.=4 max.=99999999999 statistics");
	CentroidDatafilename = "3D_ObjectCount_ResultsTable.csv";
	saveAs("Results", thresholdingOutputPath2 + CentroidDatafilename);
    close("*");

}

// LOOP --> Getting all the pixels that are surface values of the granules from 3D Object Counter with edges EXCLUDED
for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
    run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath2 = outputDir + outputFilename + "/3DSurface_Slices/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);

	slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);	
	
	run("3D Objects Counter", "threshold=128 slice=11 min.=4 max.=99999999999 exclude_objects_on_edges surfaces");
	
	run("Image Sequence... ", "dir=" + thresholdingOutputPath2 + " format=TIFF digits=2");
	resultsOutputPath = outputDir+outputFilename+"/3DSurface_Results/";
	File.makeDirectory(resultsOutputPath);
    zSliceList = getFileList(thresholdingOutputPath2);
    for(j=0;j<zSliceList.length;j++){
        open(thresholdingOutputPath2+zSliceList[j]);
        run("Image to Results"); //Results files are saved as Results0, Results1, ..
        resultsFileName = "Results"+j;
        saveAs("Results",resultsOutputPath+resultsFileName + ".csv");
        close();
    }
    close("*");
    selectWindow("Results"); 
    run("Close");
}

// LOOP --> Getting all the pixels that are surface values of the granules from 3D Object Counter with edges INCLUDED
for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
    run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath2 = outputDir + outputFilename + "/3DSurface_Slices_All/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);

	slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);		
	
	run("3D Objects Counter", "threshold=128 slice=11 min.=4 max.=99999999999 surfaces");
	
	run("Image Sequence... ", "dir=" + thresholdingOutputPath2 + " format=TIFF digits=2");
	resultsOutputPath = outputDir+outputFilename+"/3DSurface_Results_All/";
	File.makeDirectory(resultsOutputPath);
    zSliceList = getFileList(thresholdingOutputPath2);
    for(j=0;j<zSliceList.length;j++){
        open(thresholdingOutputPath2+zSliceList[j]);
        run("Image to Results"); //Results files are saved as Results0, Results1, ..
        resultsFileName = "Results"+j;
        saveAs("Results",resultsOutputPath+resultsFileName + ".csv");
        close();
    }
    close("*");
    selectWindow("Results"); 
    run("Close");
}



// LOOP --> Getting what are essentially the ERODE results for the thresholded image but using image subtraction of the granule images minus the 3D surface (edges included)
// Note: There are pixels in the first and last slice that are not counted for in the 3D surface (since they are not part of granules). 

for(i=0;i<list.length;i++){
	setBatchMode(true);
//The following actions are performed for each z-stack
    run("Import HDF5", "select=["+ inputDir+list[i] + "] datasetname=/exported_data axisorder=zyxc");
    inputFilename = substring(list[i],0,lengthOf(list[i])-22); // this is the number of characters in "Simple_Segmentation.h5" which should be at the end of the files
	outputFilename = inputFilename + "Ilastik";
	print(inputFilename);
	print(outputFilename);
	thresholdingOutputPath2 = outputDir + outputFilename + "/ErodeSubtract_Images/";
	File.makeDirectory(thresholdingOutputPath2);
	print(thresholdingOutputPath2);
	
   	imagename = inputFilename;
    masklistname = "MASK_Surface map of " + imagename;
	resultlistname = "Result of " + imagename;
	print(masklistname);
	print(resultlistname);
	
	slice_number = nSlices-2;
    run("Duplicate...", "duplicate range=1-" + slice_number);
    //run("Brightness/Contrast..."); //This next step is to change the pixel values from 1 and 2 to 0 and 255 (this is more amenable to my code later on)
	setMinAndMax(1, 2);
	call("ij.ImagePlus.setDefault16bitRange", 8);
	run("Apply LUT", "stack");
	rename(inputFilename);	
    
	run("3D Objects Counter", "threshold=128 slice=11 min.=4 max.=99999999999 surfaces");
	
	setThreshold(1, 65535, "raw");
	run("Make Binary", "method=Default background=Dark black create");
	// selectWindow(masklistname);
	imageCalculator("Subtract create stack", imagename, masklistname);
	selectWindow(resultlistname);
	
	run("Image Sequence... ", "dir=" + thresholdingOutputPath2 + " format=TIFF digits=2");
	resultsOutputPath = outputDir+outputFilename+"/ErodeSubtract_Results/";
	File.makeDirectory(resultsOutputPath);
    zSliceList = getFileList(thresholdingOutputPath2);
    for(j=0;j<zSliceList.length;j++){
        open(thresholdingOutputPath2+zSliceList[j]);
        run("Image to Results"); //Results files are saved as Results0, Results1, ..
        resultsFileName = "Results"+j;
        saveAs("Results",resultsOutputPath+resultsFileName + ".csv");
        close();
    }
    close("*");
    selectWindow("Results"); 
    run("Close");
}

setBatchMode(false);