// Automated colony count for TRA-1-60-stained wells
// by Andrew Phan

setMinAndMax(0, 160);
//originally 190

run("Apply LUT");
run("32-bit");
setAutoThreshold("MinError");

run("Convert to Mask");

run("Watershed");

run("Fill Holes");

run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=[Overlay Outlines] display include summarize record add");
close();
} ;
