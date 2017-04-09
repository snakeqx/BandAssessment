# Band Assessment Statistics Specification

## 1.  Introduction
This tool is designed for analyzing band assessment image quality for CT systems.
It try to find a way to quantify the image quality and have a statistics view.

## 2. Develop tool
The program is written in Python 3.x with several packages includes:
* pydicom -> To parse the dicom file
* matplotlib -> To draw the plot
* numpy -> To make calculation quick when dealing with image arrays
* pillow --> To generate image and do some adjustment of the image

## 3. Analyzing Method Short Introduction
The analyzing method is so-called "circular Integration".
Fake code as below:
```
DiciomImage = Image.Open("a.dcm")
Radius = DicomImage.GetRadius()

for i in range(1, Radius):
    IntegrationResult[i] = DrawCircle(DicomImage, radius).getMeanValue()
```
