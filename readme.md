# Summary

This is an application that takes two images as panorama input and stitches them together by computing a homography between the images and constructing an equirectangular image projection that stiches both images together. Overlapping pixels are blended together. For this project I built a Harris corner detector as well as a RANSAC alignment algorithm.

This project was written in C++, using the [Qt](http://qt.nokia.com/products/) framework framework. 

# Credits

For this application,[Larry Zitnick](http://research.microsoft.com/en-us/people/larryz/) at Microsoft Research wrote the majority of the UI code and the high-level architecture, and [Nat Guy](http://www.natguy.net) wrote the majority of the majority of the actual image processing and stitching code.