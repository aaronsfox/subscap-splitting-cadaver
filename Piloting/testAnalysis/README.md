# Pilot Analysis of Digitisation and Processing

This folder includes digitised images and code to process the muscle lines of action and moment arms.

## Image Used

The image selected for this pilot process was from specimen 1, 50% split, 0 degree abduction, 40 newton load and taken in the scapular plane (id: D0910_T1435).

## Digitisation Process

Mimics was used to digitise the image. The following steps were taken:

- Image directions were selected to rotate the image appropriately (i.e. top/bottom, left/right etc.) as part of the Mimics import process. ***TODO: Check orientation of these and make sure they are the same so that coordinates export in the same manner. In this case the image became 'coronal' and hence the Y coordinate was 0, and X-Z represents X-Y in 2D...***
- Points were place using Mimics 'Analyze' tab on wires representing the muscle lines of action. The naming convention for these points was firstly 'ssX' where the X relates to the line of subscapularis being referred to (i.e. 1 = top; 2 = top-middle; 3 = bottom-middle; 4 = bottom). This was followed by '_X' where the X relates to the order of the points being digitised, with incrementing numbers representing a insertion to origin ordering. For example, ss3_1 represents the digitised point closest to the insertion of the bottom-middle subscapularis line. Variable numbers of points could be digitised on each line and this is accounted for in later analysis code. In general, 2 or 3 points are fairly easy to do when the lines of action can easily be identified. These points are marked **<u>red</u>**.
- Points were also placed around the articular cartilage section of the humeral head. These points were labelled 'hhX' where the X was incrementally increased up from 1 onwards. The order of these don't necessarily matter as they will be used in a global context to fit a circle to the humeral head. These points are marked **<u>green</u>**.
- Points were similarly placed along the curved surface of the glenoid. These points were labelled 'gpX' where the X was incrementally increased up from 1 onwards. Again, the order of these don't matter as a line will be fitted to them to signify the glenoid plane. These are probably the trickiest to define and I'm less confident about the ones on this image. This approach also slightly differs to the approach of marking the most superior and inferior points of the glenoid border and conneting a line beween them to define the plane. I'm not sure which approach is better or worse. These points are marked **<u>yellow</u>**.
- A circle was fit to the circular phantom bead (labelled 'phantom'). This is not completely circle, so it's important the circular radius covers the slightly wider part so that the radius measures the 6.5mm bead diameter. This radial distance will later be used to normalise moment arm distance in the images. This circle is marked **<u>purple</u>**.

Below is a sample of the end result of the above digitisation process in Mimics.

![](sampleDigitisation.bmp)

## Export Process

From Mimics the points across muscle lines, humeral head and glenoid plane are all exported, along with the phantom bead circle. This is all done in the one text file named 'digitisationExport.txt'. To do this, the items are selected and from the right click menu select 'Export txt...'. This is then a simple process of adding the items and exporting to the appropriate folder.

The Mimics file can also be saved as a project file so it can be re-opened later for any changes.

## Data Analysis

...

...base version of Anaconda distribution with Python 3.7...