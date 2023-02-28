# Does subscapularis splitting technique during the Latarjet procedure matter?



## Data

- For each specimen, data is contained in folders labelled 'split50', 'split25_upper' and 'split25_lower' signifying the different Latarjet conditions
- Within each folder there are subfolders related to the different:
  - Joint position - abd0, abd30, abd60, abd90, ABER, APP
  - Load - 0N, 10N, 20N, 30N, 40N
  - Plane - SP, TP
- The subfolders then contain:
  - The images from the X-ray - the one used for analysis is labelled 'proc' and has had brightness/contrast adjusted for better viewing
  - The .csv files with the XY coordinates of digitised points - ss1, ss2, ss3, ss4, gp, hh, phantom, nsew, pts
    - Note that gp is only defined in the abd0 0N position for both scapular and transverse plane - and then matched to the other images
    - To cover this, there are points on the screws and ExFix pins digitised in a specific order to allow image matching (named: 'pts'). This varies per specimen, but some suggestions are:
      - For the scapular plane:
        - Top of superior screw (centre of circular portion)
        - Bottom of superior screw (most lateral aspect)
        - Top of inferior screw (centre of circular portion)
        - Bottom of inferior screw (most lateral aspect)
        - ExFix pins from left to right of image (end of pin)
      - For the transverse plane:
        - Top of the screw head (middle-top of head)
        - End of screw (middle-bottom of screw)
        - ExFix pins from left to right of image (end of pin)



## How to Digitise Images

- Open ImageJ
- Open the desired image (File > Open...)
- The image may need to be rotated (Image > Transform > Rotate)
  - For the scapular plane images, the image should be as if you are facing the shoulder from an anterior view. This means that right vs. left arm specimen will be mirror images. Whatever the case, the shaft of the humerus should be pointing down.
  - For the transverse plane images, the image should be as if you are facing the shoulder from a superior view. This means again that right vs. left arm specimen will be mirror images.
- Adjust the brightness and contrast of the image so that it is easier to view (Image > Adjust > Brightness/Contrast — or Ctrl+Shift+C)
  - The typical adjustment is to increase the Brightness and Contrast sliders
  - Hit Apply and then OK
  - This new adjusted image should be saved as a separate file (File > Save As > Tiff...). Simply save the file with **"_proc"** appended to the end
- Select the multi point digitising tool. This can be used for all the digitising selections. **The order of the points digitised across many of the following is important!** Note that in between each measurement you need to clear the points (and close the measurement box) — which is done by holding Ctrl and clicking the points.
- **North-South-East-West points**
  - Digitise four basic points at the top centre (i.e. North), bottom centre (i.e. South), right centre (i.e. East), left centre (i.e. West).
  - Use Analyze > Measure (or Ctrl+M is much easier) to record the points
  - Save this table to .csv file as **nsew.csv**
- **Phantom**
  - The circular phantom bead provides a scaling reference. We map 3 or 4 points to the outer boundary of the bead to fit a circle to. Digitise these points evenly spaced around the bead.
  - Use Analyze > Measure (or Ctrl+M is much easier) to record the points
  - Save this table to .csv file as **phantom.csv**
- **Glenoid Plane (not every image!)**
  - We digitise two points at the superior and inferior (scapular plane) / anterior and posterior points of the glenoid plane. We fit a line to these points to be representative of the glenoid in the respective plane. Digitise the two points — but only on a **good reference image** (see notes in next section on reference points). Identify and digitise the two necessary points. **You will need to make notes as to which images you have mapped the glenoid plane to!**
  - Use Analyze > Measure (or Ctrl+M is much easier) to record the points
  - Save this table to .csv file as **gp.csv**
- **Reference Image points**
  - These points are in case we cannot identify the glenoid plane on a specific image, and hence we can use these to map where the glenoid plane sits from one image to another. We digitise a series of recognisable points (e.g. screws, Ex Fix arms etc.) **of the same number and in the same order** across images.
  - Use Analyze > Measure (or Ctrl+M is much easier) to record the points
  - Save this table to .csv file as **pts.csv**
- **Subscapularis Lines of Action**
  - The four subscapularis lines need to be digitised separately and saved as their respective numbered line. In the scapular plane this is easily viewable, going from 1-4 in the superior to inferior direction. Digitise the line of action with points from the most medial attachment (i.e. at the knot) to where the line begins to curve/kink.
  - Use Analyze > Measure (or Ctrl+M is much easier) to record the points
  - Save this table to .csv file as **ssX.csv** (where X is the number 1, 2, 3 or 4)

