// This macro allows you to split the channels of all the .czi files in a directory, rename the image of each channel, and save it in a designated output directory. 
// Change the function starting at line 28 as needed for your own purposes.

//@ File (label = "Input directory", style = "directory") input
//@ File (label = "Output directory", style = "directory") output
//@ String (label = "File suffix", value = ".czi") suffix

input += File.separator;
output += File.separator;

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
list = getFileList(input);
list = Array.sort(list);
//waitForUser("");

for (i = 0; i < list.length; i++) {
	if(File.isDirectory(input + File.separator + list[i]))
		processFolder(input + File.separator + list[i]);
	if(endsWith(list[i], suffix))
		processFile(input, output, list[i]);
 }
}

// Change the suffixes of the file names as needed for your own purposes 
function processFile(input, output, file) {
setBatchMode(true);
run("Bio-Formats", "open=[" + input + file + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
run("Split Channels");
title = getTitle();
title=replace(title, "C3-", "");
title=replace(title, ".czi", "_FluoGFP.czi");
saveAs("Tiff", output + title); 
close();
title = getTitle();
title=replace(title, "C2-", "");
title=replace(title, ".czi", "_FluomApple.czi");
saveAs("Tiff", output + title); 
close();
title = getTitle();
title=replace(title, "C1-", "");
title=replace(title, ".czi", "_Fluo670.czi");
saveAs("Tiff", output + title); 
close();
// make sure to close every images befores opening the next one
run("Close All");
}

setBatchMode(false);
