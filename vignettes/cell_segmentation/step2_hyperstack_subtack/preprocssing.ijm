/*
 * Macro template to process multiple images in a folder
 */


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ String (label = "Index of selected channels", style = "string") channel
#@ String (label = "Total number of channels", style = "string") n
#@ String (label = "CD99 index", style = "integer") cd99

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
processFolder(input);
// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print(channel);
	print(input);
	print("Processing: " + input + File.separator + file);
	path = File.getParent(input + File.separator + file);
	output_file = File.getName(path);
	open(input + File.separator + file);
	selectWindow(file);
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+n+" slices=1 frames=1 display=Color");
	//tmp =  ;
	//print(tmp);
	run("Make Substack...", "channels="+channel);
	if(cd99 != -1){
		setSlice(cd99);
		run("Enhance Local Contrast (CLAHE)", "blocksize=39 histogram=256 maximum=40 mask=*None* fast_(less_accurate)");
	}
	saveAs("Tiff", output + File.separator + output_file + "_channel_substacked.tiff");
	selectWindow(file);
	saveAs("Tiff", output + File.separator + output_file + "_channel.tiff");
	close();
	close();
}
