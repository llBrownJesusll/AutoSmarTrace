
~AutoSmarTrace Manual~

For automated tracing of chains in AFM images or any bright chain-like objects on dark background with or without noise. Integrates with SmarTrace to extract chain persistence lengths.


**BEFORE USE**

Use of AutoSmarTrace requires a MATLAB distribution of 2017a or later for the Deep Learning Toolbox.

Images to be analyzed should be 512x512 in size and all be in the same directory.

Included with AutoSmarTrace are folder of sample images of simulated chains. The images for each folder contain chains generated with a persistence length specified by the folder name.


~Installation~

Move the AutoSmarTrace folder into a directory of your choice, then save the directory to the MATLAB path.

ex: |savepath('C:\Program Files\matlab\AutoSmarTrace');
    |savepath('C:\Project\Chain Tracing\AutoSmarTrace');


~Using AutoSmarTrace~

To use AutoSmarTrace, call the program with the following syntax:

|AutoSmarTrace(InputFilePath,OutputName,pixelSize,OptionalInput);


InputFilePath: string containing the path to the directory containing the images

OutputName: string containing the name for the folder created with figures of traces and name for the output data structure containing chain data.

pixelSize: size of the pixels in the image (default nanometres).

OptionalInput: optional string input; currently if this is 'EM', optimizes chain tracing for EM images, and if this is 'actin', optimizes tracing for actin filaments. Also, if it is having trouble tracing, this could be 'filter', and additional filters may be applied to improve tracing in some cases.

ex: |AutoSmarTrace('C:\Project\Image Data','Data Set 1',0.5);
ex: |AutoSmarTrace('C:\Project\Image Data','Data Set 2',3.94,'actin');


~Output~

The data structure output by AutoSmarTrace holds the spline data in the bspline structure. The points structure contains a less refined spline as well as the filepath to the image that each spline came from. In traced chain images, each chain is numbered according to its row in the output structures.
