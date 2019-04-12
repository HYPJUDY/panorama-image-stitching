# PanoramaImageStitching
Automatically stiching several individual images to generate a panorama image.
Follow the routine described in the paper [Automatic Panoramic Image Stitching using Invariant Features](http://matthewalunbrown.com/papers/ijcv2007.pdf).
For more details about theory and practice of panorama image stitching including algorithms, results and analysis, please read [my post](https://hypjudy.github.io/2017/05/10/panorama-image-stitching/). 
Any feedback is welcomed!

## Dependencies
Implement with [The CImg Library](http://cimg.eu/) and [VLFeat open source library](http://www.vlfeat.org/index.html) in C++ language.

Test on Visual Studio 2015, C++11. You can also compile in Linux with [linuxbuild.sh](https://github.com/HYPJUDY/panorama-image-stitching/blob/master/linuxbuild.sh).

Images should be of `bmp` format(much easier to convert by [ImageMagick](https://www.imagemagick.org/script/index.php)).

## Usage of VLFeat
[Download](http://www.vlfeat.org/download.html) and unpack the latest VLFeat binary distribution (I used [VLFeat 0.9.20 binary package](http://www.vlfeat.org/download/vlfeat-0.9.20-bin.tar.gz) on May 11, 2017) in a directory of your choice (e.g. `F:\vlfeat-0.9.20-bin`)
For Mac users, read this [Apple X-Code tutorial](http://www.vlfeat.org/xcode.html).
These instructions show how to setup a basic VLFeat project with Visual Studio 2015 on Windows 10 (64 bit). Instructions for other versions of Visual Studio and Windows should be similar.
1. Specify the value of the PATH environment variable. E.g. Add `F:\vlfeat-0.9.20-bin\vlfeat-0.9.20\bin\win64` to `Path`, where `win64` depends on your architecture.
2. Open VS and select `File > New > Project` and choose `General > Empty Project`.
3. Click `Project > Properties` to open Property Pages. Choose `Configuration > All Configurations`.
4. Select `Configuration Properties > C/C++ > General` and add the root path of the VLFeat folder (e.g. `F:\vlfeat-0.9.20-bin\vlfeat-0.9.20`) to `Additional Include Directories`.
5. Select `Configuration Properties > Linker > General` and add the `bin/win32` folder (even if your architecture is win64) (e.g. `F:\vlfeat-0.9.20-bin\vlfeat-0.9.20\bin\win32`) to `Additional Library Directories`. 
6. Select `Configuration Properties > Linker > Input` and add `vl.lib` to `Additional Dependencies`. 
7. Copy the `vl.dll` file in `bin/win32` folder (e.g. `F:\vlfeat-0.9.20-bin\vlfeat-0.9.20\bin\win32`) to your VS project's debug or release folder (e.g. `C:\Users\HYPJUDY\Documents\Visual Studio 2015\Projects\PanoramaImageStitching\Debug`).
8. Build and run the project (Ctrl+F5)

Reference: [Using from C](http://www.vlfeat.org/install-c.html) and [Microsoft Visual Studio tutorial](http://www.vlfeat.org/vsexpress.html).

# Results
![](https://github.com/HYPJUDY/panorama-image-stitching/blob/master/output/building-ori.jpg)
![](https://github.com/HYPJUDY/panorama-image-stitching/blob/master/output/building-cropped.bmp)

![](https://github.com/HYPJUDY/panorama-image-stitching/blob/master/output/flower-ori.jpg)
![](https://github.com/HYPJUDY/panorama-image-stitching/blob/master/output/flower-cropped.bmp)
