// Automated analysis of fluorescent images
// by Praise Sar Oo

// Getting the image file by opening and directing its location in your directory
open("a549_csrp1kd_lung5.tif");
// Selecting an image from the folder
selectImage("a549_csrp1kd_lung5.tif");

// Retrieve statistics of the image, specifically get the max pixel value
getStatistics(area, mean, min, max, std, histogram);
// Set the standard threshold value, max value of oen image after normalization
var setStandard = 39; // because of Oct 19 lung 2
lowerThreshold = 0;

if (max < setStandard){
		lowerThreshold = setStandard - max;
		run ("Add...", "value=" + lowerThreshold);
}
// Open the outline file (previously generated)
open("a549_csrp1kd_lung5.roi");
setAutoThreshold("Yen dark no-reset");
//run("Threshold...");
setThreshold(38, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
// run watershed to separate particles that are close
run("Watershed");
open("a549_csrp1kd_lung5.roi");
run("Analyze Particles...", "size=130-Infinity pixel circularity=0.3-1.00 show=Masks display exclude include summarize overlay composite");
