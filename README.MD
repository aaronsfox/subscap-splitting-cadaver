# Does subscapularis splitting technique during the Latarjet procedure matter?

*** TODO: add readme details... ***



Data folder

- For each specimen, data is contained in folders labelled 'split50', 'split25_upper' and 'split25_lower' signifying the different Latarjet conditions
- Within each folder there are subfolders related to the different:
  - Joint position - abd0, abd30, abd60, abd90, ABER, APP
  - Load - 0N, 10N, 20N, 30N, 40N
  - Plane - SP, TP
- The subfolders then contain:
  - The images from the X-ray - the one used for analysis is labelled 'proc' and has had brightness/contrast adjusted for better viewing
  - The .csv files with the XY coordinates of digitised points - ss1, ss2, ss3, ss4, gp, hh, phantom, nsew
    - Note that gp is only defined in the abd0 0N position for both scapular and transverse plane - and then matched to the other images
    - To cover this, there are points on the screws and ExFix pins digitised in a specific order to allow image matching (named: 'pts')
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

