// imagej macro for the neuron tree delineation using smc-phd
// processes given folder for .zip files (or any other extension) and apply the set of parameters  

// command example:
// java -Xmx4g -jar ~/ImageJ/ij.jar -ijpath ~/ImageJ/plugins/ -batch ~/test/run_phd.ijm ~/test/
// 
// java call, 4G ram
// ij.jar executable loading plugins from "~/ImageJ/plugins/"
// calls this very imagej macro "~/test/run_phd.ijm"  
// processes images from "~/test/"  directory (that's where the .tif or .tifs should be)

//----- parameter grid used for processing, values enetered as strings
sigmas 		= "4,6"; // preprocessing tubularity measure scales, multi scale filteiring... [2,6] value combinations commonly
						// numbers correspond to the tube diameter, so 6 would filter out the  thicker tubes, 2 would extract the thin ones
						// possible to go multiscale and take max
						// sigmas are not comma separable, a new script with sigmas = "2" or "6" would need to be run

// (*) = possible to give multiple param values separated with commas to extend the parameter grid, no spaces in between, just commas
// use them on the same preprocessed image
th 			= "0.3,0.4,0.5"; 						// (*) local maxima tolerance: real values [0,1] range, correspond to tubulairity's byte image range [0,255] 
// 0.04,0.06,0.08,0.1,0.15,0.2
											// 		the same measure is used in IJ's Find Maxima... where it is 10 by deafult 			
											//    	10/255=0.039 would be ij's default, but tweaking this value is essential
											//      0.04 +/- k * 0.01 
											// 		very small value can cause seed overflow and long computation (no mechanism to reulate that at this moment)
											//		if the data are very clean then higher values make sense as the local maxima are more distinct 
    										//      (saving midresults helps seeing how many seeds there were for each th)
no 			= "20";  						// (*) number of objects, typically, 20, can be 10, 30, maybe 40 (> objects, more processing time)
ro 			= "10";  					    // (*) number of particles resampled per object    (typically 10,12   15 slows down very much), 5 works for well suppressed background
ni 			= "10";  					    // (*) number of particles predicted per particle  (typically 10,12   15 slows down very much), 5 works for well suppressed background
step 		= "3";  						// (*) prediction step (3 is reasonable, 2 can be used as well)
kappa 		= "3"; 							// (*) von mises circular variance used for the prediction, 3 is default
ps          = "0.95"; 						// (*) 0.95 survival pty, ok to set high, around 90%
pd 			= "0.95"; 						// (*) 0.95 detectio pty, ok to set high, around 90%
krad    	= "4"; 							// (*) clustering kernel size for the mean-shift, this one is the kernel size for mean-shift, kept it as several voxels always 
kc			= "30,50"; 						// (*) Kc decay, the stopping sensitivity, so 30,50 higher the value less clutter contribution


maxiter  	= "200";						// limit iterations to <=200 iterations (this was really empirical, in case it extends too long)
maxepoch 	= "50";							// limit rounds to <=50 

savemidres  = "false";      				// "true"		: save midresults in separate folder (slows down very much)
											// "false"				: don't save the midresults
//----- done with the parameters


main_folder = getArgument;

if (main_folder=="") exit ("Need argument: folder with .tif images to process.");
if (!File.isDirectory(main_folder)) exit("Argument is not a folder!");

// add "/" to the end if it was missing
if (!endsWith(main_folder, "/")) {
	main_folder = main_folder + "/";
}

setBatchMode(true);
t_start = getTime();

files = getFileList(main_folder); // list all files and detect with parameter grid

for (i=0; i<files.length; i++ ) {
if (endsWith(files[i], ".tif")) { // select only *.tif files from the directory_path
	arg =  
	"select="+main_folder+files[i]+" "+ 
	"sigmas="+sigmas+" "+
	"th="+th+" "+
	"no="+no+" "+
	"ro="+ro+" "+
	"ni="+ni+" "+	
	"step="+step+" "+
	"kappa="+kappa+" "+
	"ps="+ps+" "+
	"pd="+pd+" "+
	"krad="+krad+" "+	
	"kc="+kc+" "+
	"maxiter="+maxiter+" "+
	"maxepoch="+maxepoch+" "+
	"savemidres=" + savemidres;

	print(arg);

	t1 = getTime();
	run("PHD", arg); // call the plugin with the parameters
	t2 = getTime();

	print("\n" + main_folder+files[i] + "\n    -> " + ((t2-t1)/60000) + " min.");

}}

t_end = getTime();
print("finished. elapsed " + ((t_end-t_start)/60000) + " min.");