// Automated colony count for NANOG-stained well scans
// by Nolan Origer
// based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7856428/

setMinAndMax(0, 1000);
run("Apply LUT");
run("Make Binary", "method=Otsu calculate black");
run("Fill Holes", "slice");
setOption("BlackBackground", true);
run("Erode", "slice");
run("Despeckle", "slice");
run("Analyze Particles...", "size=80000-Infinity circularity=0.00-1.00 show=[Overlay Masks] display summarize overlay slice");
