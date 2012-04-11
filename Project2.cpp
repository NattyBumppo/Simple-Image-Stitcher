#define _USE_MATH_DEFINES

#include "mainwindow.h"
#include "math.h"
#include "ui_mainwindow.h"
#include <QtGui>
#include "Matrix.h"
#include <time.h>

/*******************************************************************************
Draw a red cross on top of detected Harris corners.

    interestPts: array of interest points (x, y)
    numInterestPts: number of interest points
    imageDisplay: image used for drawing
*******************************************************************************/
void MainWindow::DrawInterestPoints(CIntPt *interestPts, int numInterestPts, QImage &imageDisplay)
{
   int i;
   int r, c, rd, cd;
   int imageWidth = imageDisplay.width();
   int imageHeight = imageDisplay.height();

   for(i=0;i<numInterestPts;i++)
   {
       c = (int) interestPts[i].m_X;
       r = (int) interestPts[i].m_Y;

       for(rd=-2;rd<=2;rd++)
           if(r+rd >= 0 && r+rd < imageHeight && c >= 0 && c < imageWidth)
               imageDisplay.setPixel(c, r + rd, qRgb(255, 0, 0));

       for(cd=-2;cd<=2;cd++)
           if(r >= 0 && r < imageHeight && c + cd >= 0 && c + cd < imageWidth)
               imageDisplay.setPixel(c + cd, r, qRgb(255, 0, 0));
   }
}

/*******************************************************************************
Compute simple 8d descriptors for each interest point (corner).

    image: input image
    interestPts: array of interest points
    numInterestPts: number of interest points

    If the descriptor cannot be computed, i.e. it's too close to the boundary of
    the image, its descriptor length will be set to 0.
*******************************************************************************/
void MainWindow::ComputeDescriptors(QImage image, CIntPt *interestPts, int numInterestPts)
{
    if(!numInterestPts)
	{
		return;
	}
	
	int r, c, cd, rd, i;
    int imageWidth = image.width();
    int imageHeight = image.height();
    double *buffer = new double [imageWidth * imageHeight];
    QRgb pixel;

    // Descriptor parameters
    double sigma = 2.0;
    int rad = 4;

    // Compute descriptors from green channel
    for(r=0;r<imageHeight;r++)
       for(c=0;c<imageWidth;c++)
        {
            pixel = image.pixel(c, r);
            buffer[r*imageWidth + c] = (double) qGreen(pixel);
        }

    // Blur
	SeparableGaussianBlurImage(buffer, imageWidth, imageHeight, sigma);

    // Compute the descriptor from the difference between the point sampled at its center
    // and eight points sampled around it.
    for(i=0;i<numInterestPts;i++)
    {
        int c = (int) interestPts[i].m_X;
        int r = (int) interestPts[i].m_Y;

		// If any of the sampled points falls outside of the image, reject it
        if(c >= rad && c < imageWidth - rad && r >= rad && r < imageHeight - rad)
        {
            double centerValue = buffer[(r)*imageWidth + c];
            int j = 0;

            for(rd=-1;rd<=1;rd++)
                for(cd=-1;cd<=1;cd++)
                    if(rd != 0 || cd != 0)
                {
                    interestPts[i].m_Desc[j] = buffer[(r + rd*rad)*imageWidth + c + cd*rad] - centerValue;
                    j++;
                }

            interestPts[i].m_DescSize = DESC_SIZE;
        }
        else
        {
			interestPts[i].m_DescSize = 0;
        }
    }

    delete [] buffer;
}

/*******************************************************************************
Draw a green line between matches in the two images (i.e., visualizes the potential homography).

    matches: matching points
    numMatches: number of matching points
    image1Display: image to draw matches
    image2Display: image to draw matches
*******************************************************************************/
void MainWindow::DrawMatches(CMatches *matches, int numMatches, QImage &image1Display, QImage &image2Display)
{
    int i;
    // Show matches on image
    QPainter painter;
    painter.begin(&image1Display);
    QColor green(0, 250, 0);
    QColor red(250, 0, 0);

    for(i=0;i<numMatches;i++)
    {
        painter.setPen(green);
        painter.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter.setPen(red);
        painter.drawEllipse((int) matches[i].m_X1-1, (int) matches[i].m_Y1-1, 3, 3);
    }

    QPainter painter2;
    painter2.begin(&image2Display);
    painter2.setPen(green);

    for(i=0;i<numMatches;i++)
    {
        painter2.setPen(green);
        painter2.drawLine((int) matches[i].m_X1, (int) matches[i].m_Y1, (int) matches[i].m_X2, (int) matches[i].m_Y2);
        painter2.setPen(red);
        painter2.drawEllipse((int) matches[i].m_X2-1, (int) matches[i].m_Y2-1, 3, 3);
    }

}


/*******************************************************************************
Given a set of matches, compute the "best fitting" homography.

    matches: matching points
    numMatches: number of matching points
    h: returned homography
    isForward: direction of the projection (true = image1 -> image2, false = image2 -> image1)
*******************************************************************************/
bool MainWindow::ComputeHomography(CMatches *matches, int numMatches, double h[3][3], bool isForward)
{
    int error;
    int nEq=numMatches*2;

    dmat M=newdmat(0,nEq,0,7,&error);
    dmat a=newdmat(0,7,0,0,&error);
    dmat b=newdmat(0,nEq,0,0,&error);

    double x0, y0, x1, y1;

    for (int i=0;i<nEq/2;i++)
    {
        if(isForward == false)
        {
            x0 = matches[i].m_X1;
            y0 = matches[i].m_Y1;
            x1 = matches[i].m_X2;
            y1 = matches[i].m_Y2;
        }
        else
        {
            x0 = matches[i].m_X2;
            y0 = matches[i].m_Y2;
            x1 = matches[i].m_X1;
            y1 = matches[i].m_Y1;
        }


        //Eq 1 for corrpoint
        M.el[i*2][0]=x1;
        M.el[i*2][1]=y1;
        M.el[i*2][2]=1;
        M.el[i*2][3]=0;
        M.el[i*2][4]=0;
        M.el[i*2][5]=0;
        M.el[i*2][6]=(x1*x0*-1);
        M.el[i*2][7]=(y1*x0*-1);

        b.el[i*2][0]=x0;
        //Eq 2 for corrpoint
        M.el[i*2+1][0]=0;
        M.el[i*2+1][1]=0;
        M.el[i*2+1][2]=0;
        M.el[i*2+1][3]=x1;
        M.el[i*2+1][4]=y1;
        M.el[i*2+1][5]=1;
        M.el[i*2+1][6]=(x1*y0*-1);
        M.el[i*2+1][7]=(y1*y0*-1);

        b.el[i*2+1][0]=y0;

    }
    int ret=solve_system (M,a,b);
    if (ret!=0)
    {
        freemat(M);
        freemat(a);
        freemat(b);

        return false;
    }
    else
    {
        h[0][0]= a.el[0][0];
        h[0][1]= a.el[1][0];
        h[0][2]= a.el[2][0];

        h[1][0]= a.el[3][0];
        h[1][1]= a.el[4][0];
        h[1][2]= a.el[5][0];

        h[2][0]= a.el[6][0];
        h[2][1]= a.el[7][0];
        h[2][2]= 1;
    }

    freemat(M);
    freemat(a);
    freemat(b);

    return true;
}


/*******************************************************************************
A little helper function to clamp a value between a min and a max.

	value: value to clamp
	min: range minimum
	max: range maximum
*******************************************************************************/
double clamp(double value, double min, double max)
{
	if (value < min)
	{
		value = min;
	}
	if (value > max)
	{
		value = max;
	}

	return value;
}


/*******************************************************************************
Get the L1-norm distance between two interest points.

	interestPoint1: first point
	interestPoint2: second point

	(Order doesn't matter.)

*******************************************************************************/
double getDistance(CIntPt interestPoint1, CIntPt interestPoint2)
{
	double total = 0.0;
	for(int i = 0; i < interestPoint1.m_DescSize; i++)
	{
		total += (interestPoint1.m_Desc[i] - interestPoint2.m_Desc[i]) * (interestPoint1.m_Desc[i] - interestPoint2.m_Desc[i]);
	}
	
	return sqrt(total);
}


/*******************************************************************************
Take a buffer pre-padded with empty borders of thickness "radius",
and fills these borders with zeroes.

	buffer[]: pre-padded image buffer
	radius: thickness of border in pixels
	bufferWidth: total width of buffer in pixels
	bufferHeight: total height of buffer in pixels
*******************************************************************************/
void zeroBorders(double buffer[], int radius, int bufferWidth, int bufferHeight)
{
	int r, c;
	    
	// Pad left edge with zeros
	for(r = 0; r < bufferHeight; r++)
	{
		for(c = 0; c < radius; c++)
		{
			buffer[r * bufferWidth + c] = 0.0;
		}
	}

	// Pad right edge with zeros
	for(r = 0; r < bufferHeight; r++)
	{
		for(c = bufferWidth - radius; c < bufferWidth; c++)
		{
			buffer[r * bufferWidth + c] = 0.0;
		}
	}

	// Pad top edge with zeros
	for(r = 0; r < radius; r++)
	{
		for(c = radius; c < bufferWidth - radius; c++)
		{
			buffer[r * bufferWidth + c] = 0.0;
		}
	}

	// Pad bottom edge with zeros
	for(r = bufferHeight - radius; r < bufferHeight; r++)
	{
		for(c = radius; c < bufferWidth - radius; c++)
		{
			buffer[r * bufferWidth + c] = 0.0;
		}
	}
}


/*******************************************************************************
Cull out any points that aren't local maxima in a (neighborhoodSize x neighborhoodSize) neighborhood

	image[]: input image
	imageWidth: width of image in pixels
	imageHeight: height of image in pixels
	neighborhoodSize: size of search neighborhood
*******************************************************************************/
void localMaxima(double image[], int imageWidth, int imageHeight, int neighborhoodSize)
{
	int neighborhoodRadius = (neighborhoodSize - 1) / 2;
	// Create a buffer with borders
	int bufferWidth = imageWidth + neighborhoodRadius * 2;
	int bufferHeight = imageHeight + neighborhoodRadius * 2;
	double *buffer = new double [bufferWidth * bufferHeight];
	
	// Copy image into the center of the buffer
	for(int r = 0; r < imageHeight; r++)
	{
		for(int c = 0; c < imageWidth; c++)
		{
			buffer[(r + neighborhoodRadius) * bufferWidth + (c + neighborhoodRadius)] = image[r * imageWidth + c];
		}
	}
	
	// Zero out the border
	zeroBorders(buffer, neighborhoodRadius, imageWidth, imageHeight);
	
	// For each pixel in the buffer (not counting the edges)
	for(int r = neighborhoodRadius; r < bufferHeight - neighborhoodRadius; r++)
	{
		for(int c = neighborhoodRadius; c < bufferWidth - neighborhoodRadius; c++)
		{
			// Compare to every pixel in the neighborhood (except itself)
			for(int windowRow = r - neighborhoodRadius; windowRow < r + neighborhoodRadius; windowRow++)
			{
				for(int windowCol = c - neighborhoodRadius; windowCol < c + neighborhoodRadius; windowCol++)
				{
					// If we're on our own pixel, ignore
					if((windowRow == r)&&(windowCol == c))
					{
						// Do nothing
					}
					else
					{
						// If the intensity of the middle pixel is less than the neighborhood pixel
						// to which it's being compared, make the middle pixel not count
						if(buffer[r * bufferWidth + c] <= buffer[windowRow * bufferWidth + windowCol])
						{
							image[(r - neighborhoodRadius) * imageWidth + (c - neighborhoodRadius)] = 0.0;
						}
					}
				}
			}
		}
	}

	delete [] buffer;
}

/*******************************************************************************
A test function to draw any floating-point image.

	imageToDraw[]: image to draw
	imageToDrawWidth: width of image to draw in pixels
	imageToDrawHeight: height of image to draw in pixels
	imageDisplay: display to which to draw the image

*******************************************************************************/
void testDraw(double imageToDraw[], int imageToDrawWidth, int imageToDrawHeight, QImage &imageDisplay)
{
	if(imageDisplay.width() > imageToDrawWidth)
	{
		return;
	}
	if(imageDisplay.height() > imageToDrawHeight)
	{
		return;
	}

	// Get max for scaling
	double max = 0.0;
	for(int r = 0; r < imageDisplay.height(); r++)
	{
		for(int c = 0; c < imageDisplay.width(); c++)
		{
			if((imageToDraw[r * imageToDrawWidth + c]) > max)
			{
				max = imageToDraw[r * imageToDrawWidth + c];
			}
		}
	}

	double factor = 255.0 / max;
	
	for(int r = 0; r < imageDisplay.height(); r++)
	{
		for(int c = 0; c < imageDisplay.width(); c++)
		{
			imageDisplay.setPixel(c, r, qRgb((int) imageToDraw[r * imageToDrawWidth + c] * factor, (int) imageToDraw[r * imageToDrawWidth + c] * factor, (int) imageToDraw[r * imageToDrawWidth + c] * factor));
		}
	}
}

/*******************************************************************************
Take an image and a kernel to convolve it with, and return the convolved image
(This should work for 1-D kernels as well, allowing it to be used as part of
a separable filter.) Note that this function assumes that any image passed
in has its edges pre-padded, if necessary. Padding is left intact.

	imageToConvolve[]: buffer containing image to convolve
	imageWidth: width of image in pixels
	imageHeight: height of image in pixels
	kernel[]: buffer containing convolution kernel
	kernelWidth: width of kernel in pixels
	kernelHeight: height of kernel in pixels
	verticalPadding: height of padding used (padding is necessary so that the
		convolution kernel doesn't run off the edge)
	horizontalPadding: width of padding used
	clampTo255: whether or not to clamp intensity values 0-255
	add128: whether or not to add 128 to each intensity value (in case an operation subtracted from it)

*******************************************************************************/
void convolve(double imageToConvolve[], int imageWidth, int imageHeight, double kernel[], int kernelWidth, int kernelHeight, int verticalPadding, int horizontalPadding, bool clampTo255, bool add128)
{
	int imageRow, imageCol, kernelRow, kernelCol;
	
	// Figure out the kernel radii from the dimensions.
	int vKernelRadius = (kernelHeight - 1) / 2;
	int hKernelRadius = (kernelWidth - 1) / 2;

	// Make a copy of the image to convolve for modifications
	double * imageCopy = new double [imageWidth * imageHeight];
	for(imageRow = 0; imageRow < imageHeight; imageRow++)
	{
		for(imageCol = 0; imageCol < imageWidth; imageCol++)
		{
			imageCopy[imageRow * imageWidth + imageCol] = imageToConvolve[imageRow * imageWidth + imageCol];
		}
	}

    // For each pixel in the image (except the edges)
    for(imageRow = verticalPadding; imageRow < imageHeight - verticalPadding; imageRow++)
    {
        for(imageCol = horizontalPadding; imageCol < imageWidth - horizontalPadding; imageCol++)
        {
			double newPixel = 0.0;

			// Convolve the kernel at each pixel
			for(kernelRow = 0; kernelRow < kernelHeight; kernelRow++)
			{
				for(kernelCol = 0; kernelCol < kernelWidth; kernelCol++)
				{
					double pixel;

					// Get the pixel value (adjust coordinates by horizontalPadding
					// and verticalPadding to account for edges, which aren't being used)
					pixel = imageToConvolve[(imageRow - vKernelRadius + kernelRow) * imageWidth + imageCol - hKernelRadius + kernelCol];

					// Get the value of the kernel weight
					double weight = kernel[kernelRow * kernelWidth + kernelCol];

					newPixel += weight * pixel;
				}
			}

			if(add128)
			{
				newPixel += 128.0;
			}

			if(clampTo255)
			{
				// Clamp pixel value
				newPixel = clamp(newPixel, 0.0, 255.0);
			}

			// Store convolved pixel in the image
            imageCopy[imageRow * imageWidth + imageCol] = newPixel;
		}
	}

	// Now, copy the changed copy back
	for(imageRow = 0; imageRow < imageHeight; imageRow++)
	{
		for(imageCol = 0; imageCol < imageWidth; imageCol++)
		{
			imageToConvolve[imageRow * imageWidth + imageCol] = imageCopy[imageRow * imageWidth + imageCol];
		}
	}

	delete [] imageCopy;
}


/*******************************************************************************
Blur a single channel floating point image with Gaussian weights.

    image: input and output image
    imageWidth: width of image in pixels
    imageHeight: height of image in pixels
    sigma: standard deviation of Gaussian
*******************************************************************************/
void MainWindow::SeparableGaussianBlurImage(double *image, int imageWidth, int imageHeight, double sigma)
{

	// If the sigma (standard deviation) is < 1, don't blur.
	// (Fractional sigmas result in unbalanced kernels due to
	// the Gaussian function exploding a bit.)
	if (sigma < 1)
	{
		return;
	}

	// With a kernel that using a radius of sigma * 3, 99.7% of values are accounted for.
	int radius = floor(sigma * 3 + 0.5);
	
	// Set size for convolution kernel.
	int size = radius * 2 + 1;

	// Create a buffer with edges the thickness of radius
	int bufferWidth = imageWidth + radius * 2;
	int bufferHeight = imageHeight + radius * 2;
	double *buffer = new double [bufferWidth * bufferHeight];


	// Copy image into the center of the buffer
	for(int r = 0; r < imageHeight; r++)
	{
		for(int c = 0; c < imageWidth; c++)
		{
			buffer[(r + radius) * bufferWidth + (c + radius)] = image[r * imageWidth + c];
		}
	}

	// Zero out the border
	zeroBorders(buffer, radius, bufferWidth, bufferHeight);
	
	// Compute two kernels to convolve with the image.
    double *hKernel = new double [size];
	double *vKernel = new double [size];
	double e = M_E;
	double pi = M_PI;

	// Construct kernels for horizontal and vertical blur
	for (int i = 0; i < size; i++)
	{
		// x distance is the center point (equivalent to the radius value)
		// minus the x value of the current kernel pixel
		int x = radius - i;

		// Plug and chug into the Gaussian filter 1D formula
		hKernel[i] = (1.0 / (sqrt(2.0 * pi) * sigma)) * pow(e, -((double)x * (double)x) / (2 * sigma * sigma));
		vKernel[i] = hKernel[i];
	}

	// Apply horizontal blur
	convolve(buffer, bufferWidth, bufferHeight, hKernel, size, 1, radius, radius, false, false);
	
	// Apply vertical blur
	convolve(buffer, bufferWidth, bufferHeight, vKernel, 1, size, radius, radius, false, false);

	// Copy buffer back into (normal-sized) image
	for(int r = 0; r < imageHeight; r++)
	{
		for(int c = 0; c < imageWidth; c++)
		{
			image[r * imageWidth + c] = buffer[(r + radius) * bufferWidth + (c + radius)];
		}
	}

	delete [] hKernel;
	delete [] vKernel;
	delete [] buffer;

}


/*******************************************************************************
Detect Harris corners.

    image: input image
    sigma: standard deviation of Gaussian used to blur corner detector
    thres: Threshold for detecting corners
    interestPts: returned interest points
    numInterestPts: number of interest points returned
    imageDisplay: image returned to display (for debugging)
*******************************************************************************/
void MainWindow::HarrisCornerDetector(QImage image, double sigma, double thres, CIntPt **interestPts, int &numInterestPts, QImage &imageDisplay)
{
    int r, c;
	double harrisScaleFactor = 19.0;

    int imageWidth = image.width();
    int imageHeight = image.height();
	if((!imageWidth)||(!imageHeight))
	{
		return;
	}

    double *buffer = new double [imageWidth * imageHeight];
    QRgb pixel;

    numInterestPts = 0;

    // We'll compute the corner response using just the green channel
    for(r = 0; r < imageHeight; r++)
	{
       for(c = 0; c < imageWidth; c++)
       {
            pixel = image.pixel(c, r);
			buffer[r * imageWidth + c] = (double) qGreen(pixel);
       }
	}

	
	// Define a small window for examining the image. This window will be used
	// to compute the covariance matrix and Harris response at each point.
	// Since we'll be computing a Gaussian sum of the values in the window,
	// let's size it so the radius is sigma * 3 (99.7% of values are accounted for).

	// Let's expand the buffer so that we can calculate the derivatives
	// for the original buffer data, centered in the expanded buffer.
	int expandedBufferWidth = imageWidth + 2;
	int expandedBufferHeight = imageHeight + 2;
	double *expandedBuffer = new double [expandedBufferWidth * expandedBufferHeight];
	zeroBorders(expandedBuffer, 1, expandedBufferWidth, expandedBufferHeight);

	// Now, fill in our expanded buffer with the contents of the first buffer (centered).
    for(r = 1; r < imageHeight + 1; r++)
	{
       for(c = 1; c < imageWidth + 1; c++)
       {
			expandedBuffer[r * (expandedBufferWidth) + c] = buffer[(r - 1) * imageWidth + (c - 1)];
	   }
	}
	
    // 1st-derivative kernel to convolve with the image
    double *firstDKernel = new double [3];
	firstDKernel[0] = -1.0;
	firstDKernel[1] = 0.0;
	firstDKernel[2] = 1.0;

	// Set up image copies for the x and y derivatives, as well as the arrays
	// for which we'll use them
	double * x_deriv = new double [expandedBufferWidth * expandedBufferHeight];
	double * y_deriv = new double [expandedBufferWidth * expandedBufferHeight];
	double * y_deriv_squared = new double [imageWidth * imageHeight];
	double * x_deriv_squared = new double [imageWidth * imageHeight];
	double * xy_deriv = new double [imageWidth * imageHeight];

	// Copy expanded buffer into the x- and y-derivative buffers
	for(r = 0; r < expandedBufferHeight; r++)
	{
		for(c = 0; c < expandedBufferWidth; c++)
        {
            x_deriv[r * expandedBufferWidth + c] = expandedBuffer[r * expandedBufferWidth + c];
			y_deriv[r * expandedBufferWidth + c] = expandedBuffer[r * expandedBufferWidth + c];
        }
	}

	// Calculate the x-derivative of the buffer image
	convolve(x_deriv, expandedBufferWidth, expandedBufferHeight, firstDKernel, 3, 1, 1, 1, false, false);

	// Calculate the y-derivative of the buffer image
	convolve(y_deriv, expandedBufferWidth, expandedBufferHeight, firstDKernel, 1, 3, 1, 1, false, false);

	// Calculate the x-derivative squared of the buffer image
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
        {
			x_deriv_squared[r * imageWidth + c] = x_deriv[(r + 1) * expandedBufferWidth + (c + 1)] * x_deriv[(r + 1) * expandedBufferWidth + (c + 1)];
        }
	}


	// Now blur it
	SeparableGaussianBlurImage(x_deriv_squared, imageWidth, imageHeight, sigma);

	// Calculate the y-derivative squared of the buffer image
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
        {
			y_deriv_squared[r * imageWidth + c] = y_deriv[(r + 1) * expandedBufferWidth + (c + 1)] * y_deriv[(r + 1) * expandedBufferWidth + (c + 1)];
        }
	}


	// Now blur it
	SeparableGaussianBlurImage(y_deriv_squared, imageWidth, imageHeight, sigma);

	// Calculate the x-derivative * y-derivative of the buffer image
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
        {
			xy_deriv[r * imageWidth + c] = x_deriv[(r + 1) * expandedBufferWidth + (c + 1)] * y_deriv[(r + 1) * expandedBufferWidth + (c + 1)];
        }
	}
	// Now blur it
	SeparableGaussianBlurImage(xy_deriv, imageWidth, imageHeight, sigma);

	// Set up a Harris response buffer to hold the Harris response for each pixel
	double *harrisResponseBuffer = new double [imageWidth * imageHeight];
	double dx_squared, dy_squared, dxdy, determinant, trace;

	// For each pixel in the image
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
		{
			// Get our dx^2 value
			dx_squared = x_deriv_squared[r * imageWidth + c];

			// Get our dy^2 value
			dy_squared = y_deriv_squared[r * imageWidth + c];

			// Get our dxdy value
			dxdy = xy_deriv[r * imageWidth + c];

			// Assuming a matrix like this:
			//
			//		dx^2	dxdy
			//	[					]
			//		dxdy	dy^2
			//
			
			// Find the determinant of the matrix
			determinant = dx_squared * dy_squared - dxdy * dxdy;

			// Find the trace of the matrix
			trace = dx_squared + dy_squared;

			// Finally, plug the determinant/trace for this point into the Harris response buffer
			if(trace == 0)
			{
				harrisResponseBuffer[r * imageWidth + c] = 0.0;
			}
			else
			{
				harrisResponseBuffer[r * imageWidth + c] = abs(determinant / trace) / harrisScaleFactor;
			}
		}
	}

	double min = 0.0;
	double max = 0.0;

	// Now that we've populated the Harris response buffer, let's cull out responses below the threshold
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
		{

			if((harrisResponseBuffer[r * imageWidth + c] != 0)&&(harrisResponseBuffer[r * imageWidth + c] < min))
			{
				min = harrisResponseBuffer[r * imageWidth + c];
			}
			if(harrisResponseBuffer[r * imageWidth + c] > max)
			{
				max = harrisResponseBuffer[r * imageWidth + c];
			}


			if(harrisResponseBuffer[r * imageWidth + c] < thres)
			{
				harrisResponseBuffer[r * imageWidth + c] = 0.0;
			}
		}
	}

	// Cull out any points in the Harris response that aren't local peaks in their 3x3 neighborhood
	localMaxima(harrisResponseBuffer, imageWidth, imageHeight, 5);

	// Count interest points in Harris response
	numInterestPts = 0;
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
		{
			
			if(harrisResponseBuffer[r * imageWidth + c] > 0.0)
			{
				numInterestPts++;
			}
		}
	}


	// Allocate an array in which store our interest points
	*interestPts = new CIntPt[numInterestPts];

	// Store our interest points
	int i = 0;
	for(r = 0; r < imageHeight; r++)
	{
		for(c = 0; c < imageWidth; c++)
		{
			if(harrisResponseBuffer[r * imageWidth + c] > 0.0)
			{
				(*interestPts)[i].m_X = c;
				(*interestPts)[i].m_Y = r;
				i++;
			}
		}
	}

    // Once you are done finding the interest points, display them on the image
    DrawInterestPoints(*interestPts, numInterestPts, imageDisplay);

    delete [] buffer;
	delete [] firstDKernel;
	delete [] harrisResponseBuffer;
	delete [] x_deriv;
	delete [] y_deriv;
	delete [] x_deriv_squared;
	delete [] y_deriv_squared;
	delete [] xy_deriv;
	delete [] expandedBuffer;
}


/*******************************************************************************
Find matching interest points between images.

    image1: first input image
    interestPts1: interest points corresponding to image 1
    numInterestPts1: number of interest points in image 1
    image2: second input image
    interestPts2: interest points corresponding to image 2
    numInterestPts2: number of interest points in image 2
    matches: set of matching points to be returned
    numMatches: number of matching points returned
    image1Display: image used to display matches
    image2Display: image used to display matches
*******************************************************************************/
void MainWindow::MatchInterestPoints(QImage image1, CIntPt *interestPts1, int numInterestPts1,
                             QImage image2, CIntPt *interestPts2, int numInterestPts2,
                             CMatches **matches, int &numMatches, QImage &image1Display, QImage &image2Display)
{
	
	numMatches = numInterestPts1;

    // Compute the descriptors for each interest point.
    // You can access the descriptor for each interest point using interestPts1[i].m_Desc[j].
    // If interestPts1[i].m_DescSize = 0, it was not able to compute a descriptor for that point
    ComputeDescriptors(image1, interestPts1, numInterestPts1);
    ComputeDescriptors(image2, interestPts2, numInterestPts2);

    // The position of the interest point in image 1 is (m_X1, m_Y1)
    // The position of the interest point in image 2 is (m_X2, m_Y2)
	
	*matches = new CMatches[numMatches];

	double minL2distance;
	int matchNo = 0;
	int attemptNo = 0;
	int rejectNo = 0;

	for(int image1pt = 0; image1pt < numInterestPts1; image1pt++)
	{
		// If the current point doesn't have a descriptor, skip it
		if(interestPts1[image1pt].m_DescSize == 0)
		{
			rejectNo++;
			continue;
		}

		for(int image2pt = 0; image2pt < numInterestPts2; image2pt++)
		{

			// If the current point doesn't have a descriptor, skip it
			if(interestPts2[image2pt].m_DescSize == 0)
			{
				continue;
			}

			// If it's the first match attempted, assume it's right and get the distance
			if(attemptNo == 0)
			{
				(*matches)[matchNo].m_X1 = interestPts1[image1pt].m_X;
				(*matches)[matchNo].m_Y1 = interestPts1[image1pt].m_Y;
				
				(*matches)[matchNo].m_X2 = interestPts2[image2pt].m_X;
				(*matches)[matchNo].m_Y2 = interestPts2[image2pt].m_Y;
				
				minL2distance = getDistance(interestPts1[image1pt], interestPts2[image2pt]);
				attemptNo++;
			}
			else
			{
				// If the distance is shorter, this is a better match
				if(getDistance(interestPts1[image1pt], interestPts2[image2pt]) < minL2distance)
				{
					(*matches)[matchNo].m_X1 = interestPts1[image1pt].m_X;
					(*matches)[matchNo].m_Y1 = interestPts1[image1pt].m_Y;
					
					(*matches)[matchNo].m_X2 = interestPts2[image2pt].m_X;
					(*matches)[matchNo].m_Y2 = interestPts2[image2pt].m_Y;

					minL2distance = getDistance(interestPts1[image1pt], interestPts2[image2pt]);
					attemptNo++;
				}
			}
		}
		attemptNo = 0;
		matchNo++;
	}

	numMatches = numMatches - rejectNo;

    // Draw the matches
    DrawMatches(*matches, numMatches, image1Display, image2Display);

	//delete [] (*matches);
}

/*******************************************************************************
Project a point (x1, y1) using the homography transformation h.

    (x1, y1): input point
    (x2, y2): returned point
    h: input homography used to project point
*******************************************************************************/
void Project(double x1, double y1, double &x2, double &y2, double h[3][3])
{

	// Get the product [u v w] of the following matrix multiply:
	//
	//     a  b  c        x       u
	//
	//  [  d  e  f  ] * [ y ] = [ v ]
	//
	//     g  h  1        1       w  
	//

	double u, v, w;

	u = h[0][0] * x1 + h[0][1] * y1 + h[0][2] * 1;
	v = h[1][0] * x1 + h[1][1] * y1 + h[1][2] * 1;
	w = h[2][0] * x1 + h[2][1] * y1 + h[2][2] * 1;

	x2 = u / w;
	y2 = v / w;

}

/*******************************************************************************
Count the number of inliers given a homography.  This is a helper function for RANSAC.

    h: input homography used to project points (image1 -> image2)
    matches: array of matching points
    numMatches: number of matches in the array
    inlierThreshold: maximum distance between points that are considered to be inliers

    Returns the total number of inliers.
*******************************************************************************/
int MainWindow::ComputeInlierCount(double h[3][3], CMatches *matches, int numMatches, double inlierThreshold)
{
	double projected_x, projected_y, distance;
	int numInliers = 0;

	// Project the first point in each match using the homography given,
	// and get its distance from the point it's matched to.
	for(int i = 0; i < numMatches; i++)
	{
		// Project the first point
		Project(matches[i].m_X1, matches[i].m_Y1, projected_x, projected_y, h);

		// Find the distance between the projected point and the second point
		distance = sqrt((projected_x - matches[i].m_X2) * (projected_x - matches[i].m_X2)
								+ (projected_y - matches[i].m_Y2) * (projected_y - matches[i].m_Y2));

		// If it's under the threshold, it's an inlier
		if (distance < inlierThreshold)
		{
			numInliers++;
		}
	}
    return numInliers;
}

/*******************************************************************************
Return all of the inliers, given a homography. This is a helper function for RANSAC.

    h: input homography used to project points (image1 -> image2)
    matches: array of matching points
    numMatches: number of matches in the array
    inlierThreshold: maximum distance between points that are considered to be inliers
	inliers: array of inliers calculated from the other arguments
	numInliers: size of inliers

    Returns the total number of inliers.
*******************************************************************************/
void GetInliers(double h[3][3], CMatches *matches, int numMatches, double inlierThreshold, CMatches *inliers)
{
	double projected_x, projected_y, distance;
	int numInliers = 0;

	// Project the first point in each match using the homography given,
	// and get its distance from the point it's matched to.
	for(int i = 0; i < numMatches; i++)
	{
		// Project the first point
		Project(matches[i].m_X1, matches[i].m_Y1, projected_x, projected_y, h);
		
		// Find the distance between the projected point and the second point
		distance = sqrt((projected_x - matches[i].m_X2) * (projected_x - matches[i].m_X2)
								+ (projected_y - matches[i].m_Y2) * (projected_y - matches[i].m_Y2));

		// If it's under the threshold, it's an inlier
		if (distance < inlierThreshold)
		{
			// Add inlier to list of inliers
			inliers[numInliers] = matches[i];
			// Increment inlier count
			numInliers++;
		}
	}
    

}

/*******************************************************************************
Compute homography transformation between images using RANSAC.

    matches: set of matching points between images
    numMatches: number of matching points
    numIterations: number of iterations to run RANSAC
    inlierThreshold: maximum distance between points that are considered to be inliers
    hom: returned homography transformation (image1 -> image2)
    homInv: returned inverse homography transformation (image2 -> image1)
    image1Display: first image used to display matches
    image2Display: second image used to display matches
*******************************************************************************/
void MainWindow::RANSAC(CMatches *matches, int numMatches, int numIterations, double inlierThreshold,
                        double hom[3][3], double homInv[3][3], QImage &image1Display, QImage &image2Display)
{
	// We'll be comparing groups of 4 points
	#define MATCH_GROUP_SIZE 4
	
	// If there are fewer than matchGroupSize matches, this won't work, so return.
	if(numMatches < MATCH_GROUP_SIZE)
	{
		return;
	}
	
	CMatches potentialInliers[MATCH_GROUP_SIZE];
	int numInliers, maxInliers, randomMatchID;
	int usedMatchIDs[MATCH_GROUP_SIZE];
	double potentialInlierHom[3][3];
	double bestHom[3][3];
	
	for(int i = 0; i < MATCH_GROUP_SIZE; i++)
	{
		usedMatchIDs[i] = -1;
	}

	// Initialize random seed
	srand( time(NULL) );	

	maxInliers = -1;

	for (int iter = 0; iter < numIterations; iter++)
	{
		// Randomly select 4 pairs of potentially matching points from matches
		for(int i = 0; i < MATCH_GROUP_SIZE; i++)
		{
			// Make sure that the match selected is unique
			bool matchUnique = false;
			while(!matchUnique)
			{
				matchUnique = true;
				randomMatchID = rand() % numMatches;
			
				// Compare generated match to every other match,
				// and make sure we're not reusing one
				for(int j = 0; j < i; j++)
				{
					
					// Check for a match with a previously used match
					if(randomMatchID == usedMatchIDs[j])
					{
						matchUnique = false;
						continue;
					}
				}
			}
			
			// Once we've reached this point, we know that we have a unique match,
			// so let's add the randomly generated match to the array of potential inliers
			potentialInliers[i] = matches[randomMatchID];
			
			// Make sure we don't use the same match again
			usedMatchIDs[i] = randomMatchID;
		}

		// Now that we have four unique, randomly selected matches,
		// compute the homography relating the four selected matches
		ComputeHomography(potentialInliers, MATCH_GROUP_SIZE, potentialInlierHom, true);

		// Using the computed homography, compute the number of inliers against all of the matches
		numInliers = ComputeInlierCount(potentialInlierHom, matches, numMatches, inlierThreshold);

		// If this homography produces the highest number of inliers, store it as the best homography
		if(numInliers > maxInliers)
		{
			maxInliers = numInliers;
			for(int i = 0; i < 3; i++)
			{
				for(int j = 0; j < 3; j++)
				{
					bestHom[i][j] = potentialInlierHom[i][j];
				}
			}
		}
	}

	CMatches *inliers = new CMatches[maxInliers];

    // Given the highest scoring homography, once again find all the inliers
	GetInliers(bestHom, matches, numMatches, inlierThreshold, inliers);
	numInliers = ComputeInlierCount(bestHom, matches, numMatches, inlierThreshold);

	// Compute a new refined homography using all of the inliers
	ComputeHomography(inliers, maxInliers, hom, true);

	// Compute an inverse homography as well
	ComputeHomography(inliers, maxInliers, homInv, false);
	
	// Display the inlier matches
    DrawMatches(inliers, numInliers, image1Display, image2Display);

	delete inliers;
}

/*******************************************************************************
Bilinearly interpolate image (helper function for Stitch)

    image: input image
    (x, y): location to interpolate
    rgb: returned color values
*******************************************************************************/
bool MainWindow::BilinearInterpolation(QImage *image, double newX, double newY, double rgb[3])
{
	// Computing values for pixels using this setup:
	//
	//    f(x,y)              f(x+1, y)
	//      p1-------------------p2
	//      |      f(x+a, y+b)    |
	//      |            O        |
	//      |                     |
	//      |                     |
	//      p3-------------------p4
	//  f(x,y+1)              f(x+1,y+1)
	//

	int x = floor(newX);
	int y = floor(newY);
	double a = newX - x;
	double b = newY - y;

	// Get colors for pixels in the diagram above
	QRgb p1 = image->pixel(x, y);
	QRgb p2 = image->pixel(x + 1, y);
	QRgb p3 = image->pixel(x, y + 1);
	QRgb p4 = image->pixel(x + 1, y + 1);

	// Use interpolation formula to calculate the color channels for f(x+a, y+b)
	rgb[0] = (1.0 - a) * (1.0 - b) * (double) qRed(p1)
		+ a * (1 - b) * (double) qRed(p2)
		+ (1 - a) * b * (double) qRed(p3)
		+ a * b * (double) qRed(p4);
	rgb[1] = (1.0 - a) * (1.0 - b) * (double) qGreen(p1)
		+ a * (1 - b) * (double) qGreen(p2)
		+ (1 - a) * b * (double) qGreen(p3)
		+ a * b * (double) qGreen(p4);
	rgb[2] = (1.0 - a) * (1.0 - b) * (double) qBlue(p1)
		+ a * (1 - b) * (double) qBlue(p2)
		+ (1 - a) * b * (double) qBlue(p3)
		+ a * b * (double) qBlue(p4);

    return true;
}


/*******************************************************************************
Stitch together two images using the homography transformation

    image1: first input image
    image2: second input image
    hom: homography transformation (image1 -> image2)
    homInv: inverse homography transformation (image2 -> image1)
    stitchedImage: returned stitched image
*******************************************************************************/
void MainWindow::Stitch(QImage image1, QImage image2, double hom[3][3], double homInv[3][3], QImage &stitchedImage)
{
    // Width and height of stitchedImage
    int sWidth, sHeight;

	// To compute the width and height of stitchedImage, let's first get the corners of image2
	
	// Corners are like so (image2):
	//
	//    A           B
	//
	//
	//
	//    C           D
	//
	// (x, y) is for the original coordinates in image2's corners, while
	// (x_proj, y_proj) is the projected coordinates of image2's corners projected into image1

	double image2CornerA_x = 0.0;
	double image2CornerA_y = 0.0;
	double image2CornerB_x = image2.width();
	double image2CornerB_y = 0.0;
	double image2CornerC_x = 0.0;
	double image2CornerC_y = image2.height();
	double image2CornerD_x = image2.width();
	double image2CornerD_y = image2.height();

	double image2CornerA_x_proj, image2CornerA_y_proj,
		image2CornerB_x_proj, image2CornerB_y_proj,
		image2CornerC_x_proj, image2CornerC_y_proj,
		image2CornerD_x_proj, image2CornerD_y_proj;

	// Project points into image1
	Project(image2CornerA_x, image2CornerA_y, image2CornerA_x_proj, image2CornerA_y_proj, homInv);
	Project(image2CornerB_x, image2CornerB_y, image2CornerB_x_proj, image2CornerB_y_proj, homInv);
	Project(image2CornerC_x, image2CornerC_y, image2CornerC_x_proj, image2CornerC_y_proj, homInv);
	Project(image2CornerD_x, image2CornerD_y, image2CornerD_x_proj, image2CornerD_y_proj, homInv);

	// See how large we'll need to make the combined image
	sWidth = (int) ceil(max((double) image1.width(), max(image2CornerB_x_proj, image2CornerD_x_proj))
		- min(0.0, min(image2CornerA_x_proj, image2CornerC_x_proj)));
	sHeight = (int) ceil(max((double) image1.height(), max(image2CornerC_y_proj, image2CornerD_y_proj))
		- min(0.0, min(image2CornerA_y_proj, image2CornerB_y_proj)));

    stitchedImage = QImage(sWidth, sHeight, QImage::Format_RGB32);
    stitchedImage.fill(qRgb(0,0,0));

    // Copy image 1 into stitched image at offset determined by difference of dimensions
	// between it and the stitched image if image1 is to the right/below image2;
	// if image1 is to the left/above image2, then place it at 0 for the x-/y-axis
	
	int x_offset, y_offset;

	// If image2's farthest left corner is to the left of image1, give image1 an x offset
	if(min(image2CornerA_x_proj, image2CornerC_x_proj) < 0.0)
	{
		x_offset = stitchedImage.width() - image1.width();
	}
	else
	{
		x_offset = 0;
	}

	// If image2's farthest top corner is above image1, give image1 a y offset
	if(min(image2CornerA_y_proj, image2CornerB_y_proj) < 0.0)
	{
		y_offset = stitchedImage.height() - image1.height();
	}
	else
	{
		y_offset = 0;
	}
	
	for(int row = y_offset; row < image1.height() + y_offset; row++)
	{
		for(int col = x_offset; col < image1.width() + x_offset; col++)
		{
			stitchedImage.setPixel(col, row, image1.pixel(col - x_offset, row - y_offset));
		}
	}

	double stitchedImage_x_proj, stitchedImage_y_proj;

	// Place image 2 into stitched image
	for(int row = 0; row < stitchedImage.height(); row++)
	{
		for(int col = 0; col < stitchedImage.width(); col++)
		{
			// Project this point onto image2
			Project(col - x_offset, row - y_offset, stitchedImage_x_proj, stitchedImage_y_proj, hom);
			
			// If the projected point is within image2's boundaries,
			// add the pixel value to the stitched image
			if( (0.0 < stitchedImage_x_proj) &&
				(stitchedImage_x_proj < image2.width()) &&
				(0.0 < stitchedImage_y_proj) &&
				(stitchedImage_y_proj < image2.height()))
			{
				double color[3];
				BilinearInterpolation(&image2, stitchedImage_x_proj, stitchedImage_y_proj, color);
				// Modify colors to blend edges, if image1 is already drawn here
				if(stitchedImage.pixel(col, row) != qRgb(0, 0, 0))
				{
					// Weight each image equally
					color[0] = 0.5 * color[0] + 0.5 * qRed(stitchedImage.pixel(col, row));
					color[1] = 0.5 * color[1] + 0.5 * qGreen(stitchedImage.pixel(col, row));
					color[2] = 0.5 * color[2] + 0.5 * qBlue(stitchedImage.pixel(col, row));
				}
				else
				{
					// Just use the image2 colors as-is
				}
				stitchedImage.setPixel(col, row, qRgb((int) color[0], (int) color[1], (int) color[2]));
			}
		}
	}
}