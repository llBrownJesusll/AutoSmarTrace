
~AutoSmarTrace Manual~

For automated tracing of chains in AFM images or any bright chain-like objects on dark background with or without noise.


**BEFORE USE**

Use of AutoSmarTrace requires a MATLAB distribution of 2017a or later for the Deep Learning Toolbox.

Images to be analyzed should be 512x512 in size and all be in the same directory.


~Installation~

Extract the AutoSmarTrace.zip file into a directory of your choice, then save the directory to the MATLAB path.

ex: savepath('C:\Program Files\matlab\AutoSmarTrace');
    savepath('C:\Project\Chain Tracing\AutoSmarTrace');


~Using AutoSmarTrace~

To use AutoSmarTrace, call the program with the following syntax:

|AutoSmarTrace(InputFilePath,OutputName,pixelSize,OptionalInput);


InputFilePath: string containing the path to the directory containing the images

OutputName: string containing the name for the folder created with figures of traces and name for the output data structure containing chain data.

pixelSize: size of the pixels in the image (default nanometres).

OptionalInput: optional string input; currently if this is 'EM', optimizes chain tracing for EM images.

ex: AutoSmarTrace('C:\Project\Image Data','Data Set 1',2.0);


~Output~

The data structure output by AutoSmarTrace holds the spline data in the bspline structure. The points structure contains a less refined spline as well as the filepath to the image that each spline came from. In traced chain images, each chain is numbered according to its row in the output structures.