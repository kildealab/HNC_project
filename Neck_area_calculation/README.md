# HNC project

This code was written to calculate the area/volume of a region of the neck using DICOM files from a Radiotherapy treatment. This code was created as part of Aixa Andrade's Master thesis project to graduate from McGill University: "Exploration of potential parameters that influence when to replan Head and Neck Cancer patients during Radiotherapy". Hossein Naseri helped Aixa Andrade with the implementation of 3D Registration between CBCT and CT. The supervisor of Aixa Andrade is John Kildea.  

The input files are DICOM images stored in nested folders, organized by ID. Each ID folder should have the following information (this can be modified according to your own purposes but the code is structureD to read certain format): 

1. CBCT folders named as CBCT_ID_Plan_Replan_FractionNumberPlan_FractionNumberCourse. In each CBCT folder, there must be a DICOM Registration file (RE). 
2. A folder with each RP and RS DICOM file.

The output is a dataframe with area/volume of each 3D DICOM image. The area is the cross section area of one CBCT slice. The volume corresponds to the volume of 5 CBCT slices. The neck area works with 3D registration (CBCT-CT), while the volume only works with alineation in the z axis. (It is just a prototype)
The slices of interest are located around the centre of the submandibular gland in the z direction.

The code also returns plots of each image processed. Before you run the code, you need to create the "plots_folder" that store the images.

This project contains one module of functions funcs.py, that are imported in the neck_area.py and neck_volume.py. 


To succesfully use the neck_area.py software you need to define the "images_path", the "plots_folder_path" and the "patients" file. The neck_volume.py requires the "images_path" and the ID of the folder to be examined, while funcs.py only requires the definition of the "plots_folder_path".

The "images path" is the general path that contains all the nested folders. The "plots_folder_path" can be created arbitrarily. The "patients" file needs to contain the IDs of each radiotherapy treatment. 

A paths file must contain the paths to the mentioned folders.



References:

Chloe, Litrico; John, Kildea; Haley, Patrick. 2018. Radiomics for Prostate Cancer: https://www.authorea.com/users/232204/articles/296788-radiomics-for-prostate-cancer


DICOM INNOLITICS: https://dicom.innolitics.com/ciods

DICOM STANDARD: https://www.dicomstandard.org/current

OpenCV: https://opencv.org/about/

OS: https://docs.python.org/3/library/os.html

PyDICOM: https://pydicom.github.io/

SCIPY: https://www.scipy.org/docs.html


Stackoverflow: https://stackoverflow.com/questions
